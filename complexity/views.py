import base64
from io import BytesIO
from django.shortcuts import render

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .utils import (
    validate_sequence,
    calculate_linguistic_complexity,
    calculate_entropy_profile,
    calculate_monomer_profile,
    calculate_gc_profile,
    get_statistics,
    parse_bed_file,
    find_local_minima,
    get_sequence_region,
    find_short_tandem_repeats,
    count_ling_complexity
)


def index(request):
    """Главная страница + обработка формы"""
    
    if request.method == 'POST':
        
        if 'fasta_file' not in request.FILES:
            return render(request, 'complexity/index.html', {'error': 'Пожалуйста, выберите FASTA файл'})
        
        fasta_file = request.FILES['fasta_file']
        method = request.POST.get('method', 'linguistic')
        window_size = int(request.POST.get('window_size', 50))
        step = int(request.POST.get('step', 100))
        alphabet_size = int(request.POST.get('alphabet_size', 4))
        bed_file = request.FILES.get('bed_file')
        
        try:
            min_repeat_len = int(request.POST.get('min_repeat_len', 4))
        except ValueError:
            min_repeat_len = 4
        
        try:
            # Чтение FASTA файла
            content = fasta_file.read().decode('utf-8')
            lines = content.split('\n')
            sequence = ''
            header = ''
            for line in lines:
                if line.startswith('>'):
                    if not header:
                        header = line[1:].strip()
                else:
                    sequence += line.strip().upper()
            
            if not sequence:
                return render(request, 'complexity/index.html', {'error': 'Файл пуст или имеет неверный формат'})
            
            is_valid, ambiguous_found, msg = validate_sequence(sequence)
            if not is_valid:
                return render(request, 'complexity/index.html', {'error': msg})
            
            if len(sequence) < window_size:
                return render(request, 'complexity/index.html', {'error': f'Длина последовательности ({len(sequence)}) меньше длины окна ({window_size})'})
            
            # Названия методов
            method_names = {
                'linguistic': 'Лингвистическая сложность',
                'shannon': 'Энтропия Шеннона',
                'shannon2': 'Энтропия Шеннона (2-й порядок)',
                'monomer': 'Мера мономеров',
                'gc': 'GC-содержание'
            }
            
            method_name = method_names.get(method, method)
            
            # Расчёт профиля
            profile = None
            threshold = 0.7
            ylabel = ''
            color = 'blue'
            
            if method == 'linguistic':
                profile = calculate_linguistic_complexity(sequence, window_size, step, alphabet_size)
                ylabel = 'Лингвистическая сложность'
                color = 'blue'
            elif method == 'shannon':
                profile = calculate_entropy_profile(sequence, window_size, step, 1, alphabet_size)
                ylabel = 'Энтропия Шеннона'
                color = 'green'
            elif method == 'shannon2':
                profile = calculate_entropy_profile(sequence, window_size, step, 2, alphabet_size)
                ylabel = 'Энтропия Шеннона (2-й порядок)'
                color = 'teal'
            elif method == 'monomer':
                profile = calculate_monomer_profile(sequence, window_size, step)
                ylabel = 'Мера мономеров'
                color = 'red'
            elif method == 'gc':
                profile = calculate_gc_profile(sequence, window_size, step)
                ylabel = 'GC-содержание'
                color = 'orange'
                threshold = None
            else:
                return render(request, 'complexity/index.html', {'error': 'Неизвестный метод анализа'})
            
            if not profile:
                return render(request, 'complexity/index.html', {'error': 'Ошибка при расчете профиля'})
            
            # Парсинг BED
            genes = []
            if bed_file:
                try:
                    genes = parse_bed_file(bed_file)
                    print(f"Загружено генов из BED: {len(genes)}")
                except Exception as e:
                    print(f"Ошибка при парсинге BED: {e}")
            
            # Статистика
            stats = get_statistics(profile)
            
            # Поиск участков низкой сложности
            all_low_points = [(pos, val) for pos, val in profile if threshold is not None and val < threshold]
            
            # Поиск локальных минимумов
            local_minima = []
            enhanced_minima = []
            if threshold is not None:
                local_minima = find_local_minima(profile, window=5, threshold=threshold)
                for pos, val, idx in local_minima:
                    region, start, end = get_sequence_region(sequence, pos, window=30)
                    enhanced_minima.append({
                        'position': pos,
                        'value': val,
                        'sequence_region': region,
                        'region_start': start,
                        'region_end': end
                    })
            
            # ========== ПОИСК ТАНДЕМНЫХ ПОВТОРОВ ==========
            all_repeats = find_short_tandem_repeats(sequence, min_len=2)
            
            significant_repeats = []
            for r in all_repeats:
                total_len = r['end'] - r['start']
                if total_len >= min_repeat_len:
                    # Проверяем, находится ли повтор в гене (только если есть BED)
                    in_gene = False
                    if genes:
                        for gene in genes:
                            if gene['start'] <= r['start'] <= gene['end']:
                                in_gene = True
                                break
                    
                    # Вычисляем сложность окна вокруг повтора
                    context_start = max(0, r['start'] - 30)
                    context_end = min(len(sequence), r['end'] + 30)
                    context = sequence[context_start:context_end]
                    if len(context) >= 10:
                        r['complexity'] = round(count_ling_complexity(context, 4), 4)
                    else:
                        r['complexity'] = None
                    
                    r['total_len'] = total_len
                    r['in_gene'] = in_gene
                    significant_repeats.append(r)
            
            # Вывод в консоль для отладки
            print("=" * 60)
            print(f"🔍 ПОИСК ТАНДЕМНЫХ ПОВТОРОВ")
            print(f"Длина генома: {len(sequence)} нт")
            print(f"Всего найдено повторов: {len(all_repeats)}")
            print(f"Значимых повторов (>= {min_repeat_len} нт): {len(significant_repeats)}")
            if significant_repeats:
                print("Первые 10 повторов:")
                for r in significant_repeats[:10]:
                    print(f"  {r['start']}-{r['end']}: {r['type']} {r['motif']} (длина {r['total_len']} нт, в гене: {r['in_gene']})")
            print("=" * 60)
            
            # Находим глобальный минимум
            global_min_pos = None
            global_min_val = 1.0
            for pos, val in profile:
                if val < global_min_val:
                    global_min_val = val
                    global_min_pos = pos
            
            # Находим ген для глобального минимума
            global_min_gene = None
            for gene in genes:
                if gene['start'] <= global_min_pos <= gene['end']:
                    global_min_gene = gene
                    break
            
            # Генерация графика
            positions = [p[0] for p in profile]
            values = [v[1] for v in profile]
            
            fig, ax = plt.subplots(figsize=(14, 6))
            ax.plot(positions, values, '-', linewidth=1, color=color, label=method_name)
            
            if threshold is not None:
                ax.axhline(y=threshold, color='r', linestyle='--', label=f'Порог ({threshold})')
                ax.fill_between(positions, values, threshold, 
                               where=(np.array(values) < threshold), 
                               color='red', alpha=0.3)
                
                for pos, val, _ in local_minima[:10]:
                    ax.plot(pos, val, 'ro', markersize=5)
            
            if global_min_pos:
                ax.plot(global_min_pos, global_min_val, 'r*', markersize=12, 
                       label=f'Глобальный минимум: {global_min_val:.3f}')
                ax.annotate(f'Минимум\n{global_min_val:.3f}', 
                           xy=(global_min_pos, global_min_val), 
                           xytext=(10, -20), textcoords='offset points',
                           fontsize=9, bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.8))
            
            # Отметка генов из BED
            if genes:
                y_min, y_max = ax.get_ylim()
                for gene in genes:
                    if gene['start'] < len(sequence) and gene['end'] < len(sequence):
                        center = (gene['start'] + gene['end']) / 2
                        ax.annotate(gene['name'], xy=(center, y_max * 0.95), 
                                   ha='center', fontsize=8, rotation=45,
                                   bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.7))
                        ax.axvspan(gene['start'], gene['end'], alpha=0.2, color='green')
            
            ax.set_xlabel('Позиция в геноме (нт)')
            ax.set_ylabel(ylabel)
            title_text = f'{method_name}'
            if header:
                title_text += f' | {header[:80]}...' if len(header) > 80 else f' | {header}'
            ax.set_title(title_text)
            ax.legend(loc='upper right')
            ax.grid(True, alpha=0.3)
            if method != 'gc':
                ax.set_ylim(0, 1.05)
            
            plt.tight_layout()
            
            buffer = BytesIO()
            plt.savefig(buffer, format='png', dpi=100)
            buffer.seek(0)
            image_base64 = base64.b64encode(buffer.getvalue()).decode()
            buffer.close()
            plt.close()
            
            # Биологическая интерпретация
            biological_interpretation = ""
            if threshold is not None and global_min_val < threshold:
                if global_min_gene:
                    biological_interpretation = f"""
                    <strong>Биологическая интерпретация:</strong><br>
                    Найден участок низкой сложности (значение {global_min_val:.4f}) в гене <strong>{global_min_gene['name']}</strong> 
                    на позиции {global_min_pos} нт.
                    """
                else:
                    biological_interpretation = f"""
                    <strong>Биологическая интерпретация:</strong><br>
                    Найден участок низкой сложности (значение {global_min_val:.4f}) на позиции {global_min_pos} нт.
                    """
            
            context = {
                'show_result': True,
                'image_base64': image_base64,
                'header': header[:200] if header else 'Последовательность',
                'method_name': method_name,
                'avg_value': round(stats['avg'], 4),
                'min_value': round(stats['min'], 4),
                'max_value': round(stats['max'], 4),
                'low_regions_count': len(all_low_points),
                'low_regions': all_low_points[:30],
                'local_minima_count': len(local_minima),
                'enhanced_minima': enhanced_minima[:15],
                'short_repeats': significant_repeats,
                'global_min_pos': global_min_pos,
                'global_min_val': round(global_min_val, 4),
                'global_min_gene': global_min_gene,
                'biological_interpretation': biological_interpretation,
                'sequence_length': len(sequence),
                'window_size': window_size,
                'step': step,
                'genes': genes,
                'method': method,
                'threshold': threshold,
                'min_repeat_len': min_repeat_len,
            }
            
            return render(request, 'complexity/index.html', context)
            
        except Exception as e:
            import traceback
            traceback.print_exc()
            return render(request, 'complexity/index.html', {'error': f'Ошибка: {str(e)}'})
    
    return render(request, 'complexity/index.html')