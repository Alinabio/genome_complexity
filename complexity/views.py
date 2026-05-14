import base64
import json
import csv
from io import BytesIO
from django.shortcuts import render, redirect
from django.http import HttpResponse

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
            
            request.session['profile_data'] = profile
            request.session['method_name'] = method_name
            request.session['window_size'] = window_size
            request.session['step'] = step
            
            # Парсинг BED
            genes = []
            if bed_file:
                try:
                    genes = parse_bed_file(bed_file)
                    print(f"Загружено генов из BED: {len(genes)}")
                except Exception as e:
                    print(f"Ошибка при парсинге BED: {e}")
            
            stats = get_statistics(profile)
            
            all_low_points = [(pos, val) for pos, val in profile if threshold is not None and val < threshold]
            
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
                    in_gene = False
                    if genes:
                        for gene in genes:
                            if gene['start'] <= r['start'] <= gene['end']:
                                in_gene = True
                                break
                    
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
            
            print("=" * 60)
            print(f"🔍 ПОИСК ТАНДЕМНЫХ ПОВТОРОВ")
            print(f"Длина генома: {len(sequence)} нт")
            print(f"Всего найдено повторов: {len(all_repeats)}")
            print(f"Значимых повторов (>= {min_repeat_len} нт): {len(significant_repeats)}")
            print("=" * 60)
            
            # Глобальный минимум (исключая поли-A хвост)
            global_min_pos = None
            global_min_val = 1.0
            for pos, val in profile:
                if val < global_min_val and pos < len(sequence) - 500:
                    global_min_val = val
                    global_min_pos = pos
            
            global_min_gene = None
            for gene in genes:
                if gene['start'] <= global_min_pos <= gene['end']:
                    global_min_gene = gene
                    break
            
            # ========== ГЕНЕРАЦИЯ ГРАФИКА ==========
            positions = [p[0] for p in profile]
            values = [v[1] for v in profile]
            
            # Создаём фигуру с тремя областями: график, ось X, гены
            fig = plt.figure(figsize=(14, 9))
            
            # Основная ось для графика сложности
            ax_main = plt.axes([0.08, 0.30, 0.88, 0.60])
            
            # Ось для генов (под основной, ниже подписи оси X)
            ax_genes = plt.axes([0.08, 0.10, 0.88, 0.12])
            
            # ===== ОСНОВНОЙ ГРАФИК =====
            ax_main.plot(positions, values, '-', linewidth=1.5, color=color, label=method_name, zorder=2)
            
            if threshold is not None:
                ax_main.axhline(y=threshold, color='r', linestyle='--', linewidth=1, label=f'Порог ({threshold})', zorder=1)
                ax_main.fill_between(positions, values, threshold, 
                                    where=(np.array(values) < threshold), 
                                    color='red', alpha=0.25, zorder=1)
            
            if global_min_pos:
                ax_main.plot(global_min_pos, global_min_val, 'r*', markersize=12, zorder=3)
            
            ax_main.set_ylabel(ylabel, fontsize=11)
            title_text = f'{method_name}'
            if header:
                title_text += f' | {header[:80]}...' if len(header) > 80 else f' | {header}'
            ax_main.set_title(title_text, fontsize=12)
            ax_main.grid(True, alpha=0.2, linestyle='--')
            ax_main.set_xlim(min(positions), max(positions))
            
            if method != 'gc':
                ax_main.set_ylim(bottom=-0.05, top=1.05)
            
            # Подпись оси X — на основной оси, но гены будут ниже
            ax_main.set_xlabel('Позиция в геноме (нт)', fontsize=11)
            
            # ===== ОТДЕЛЬНАЯ ПАНЕЛЬ ДЛЯ ГЕНОВ (БЕЗ ОСИ X) =====
            ax_genes.set_xlim(ax_main.get_xlim())
            ax_genes.set_ylim(0, 1)
            ax_genes.set_yticks([])
            ax_genes.set_xticks([])
            ax_genes.set_xlabel('')
            ax_genes.set_ylabel('Гены', fontsize=8, rotation=0, ha='right', va='center')
            ax_genes.yaxis.set_label_coords(-0.05, 0.5)
            
            for spine in ax_genes.spines.values():
                spine.set_visible(False)
            
            if genes:
                genes_sorted = sorted(genes, key=lambda x: x['start'])
                bar_height = 0.3
                y_pos = 0.5
                
                for i, gene in enumerate(genes_sorted):
                    if gene['start'] < len(sequence) and gene['end'] < len(sequence):
                        start = max(ax_main.get_xlim()[0], gene['start'])
                        end = min(ax_main.get_xlim()[1], gene['end'])
                        color_gene = '#2ecc71' if i % 2 == 0 else '#27ae60'
                        
                        rect = plt.Rectangle((start, y_pos - bar_height/2), 
                                            end - start, bar_height,
                                            facecolor=color_gene, alpha=0.7, edgecolor='darkgreen', linewidth=0.5)
                        ax_genes.add_patch(rect)
                        
                        center = (start + end) / 2
                        # Подписи генов — вертикальные
                        ax_genes.text(center, y_pos + bar_height/2 + 0.05, gene['name'], 
                                    ha='center', va='bottom', fontsize=6, 
                                    rotation=90, fontweight='bold')
            else:
                ax_genes.text(ax_main.get_xlim()[0] + (ax_main.get_xlim()[1] - ax_main.get_xlim()[0])/2, 0.5, 
                            'BED-файл не загружен — аннотация генов отсутствует',
                            ha='center', va='center', fontsize=8, style='italic', color='gray')
            
            # Информация о глобальном минимуме — внизу
            info_text = f"📊 {method_name} | Порог: {threshold} | "
            info_text += f"Глобальный минимум: {global_min_val:.3f} на позиции {global_min_pos} нт"
            if global_min_gene:
                info_text += f" (ген {global_min_gene['name']})"
            
            fig.text(0.5, 0.03, info_text, ha='center', va='bottom', 
                    fontsize=9, style='italic', 
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="#f0f0f0", alpha=0.8))
            
            plt.tight_layout()
            plt.subplots_adjust(bottom=0.12)
            
            buffer = BytesIO()
            plt.savefig(buffer, format='png', dpi=100)
            buffer.seek(0)
            image_base64 = base64.b64encode(buffer.getvalue()).decode()
            buffer.close()
            plt.close()
            
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
                'method': method,
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
                'threshold': threshold,
                'min_repeat_len': min_repeat_len,
                'profile_data_json': json.dumps(profile),
            }
            
            return render(request, 'complexity/index.html', context)
            
        except Exception as e:
            import traceback
            traceback.print_exc()
            return render(request, 'complexity/index.html', {'error': f'Ошибка: {str(e)}'})
    
    return render(request, 'complexity/index.html')


def export_results(request):
    """Экспорт результатов анализа в CSV файл"""
    
    if request.method == 'POST':
        method = request.POST.get('method', 'linguistic')
        window_size = request.POST.get('window_size', '50')
        step = request.POST.get('step', '100')
        
        profile_json = request.POST.get('profile_data', '')
        profile_data = []
        
        if profile_json:
            try:
                profile_data = json.loads(profile_json)
            except:
                profile_data = []
        
        if not profile_data:
            profile_data = request.session.get('profile_data', [])
        
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = f'attachment; filename="complexity_results_{method}.csv"'
        
        writer = csv.writer(response)
        writer.writerow(['# Результаты анализа лингвистической сложности'])
        writer.writerow([f'# Метод: {method}'])
        writer.writerow([f'# Длина окна: {window_size}'])
        writer.writerow([f'# Шаг: {step}'])
        writer.writerow([f'# Количество точек: {len(profile_data)}'])
        writer.writerow([])
        writer.writerow(['Позиция (центр окна)', 'Значение сложности'])
        
        for item in profile_data:
            if isinstance(item, (list, tuple)) and len(item) >= 2:
                writer.writerow([item[0], item[1]])
            elif isinstance(item, dict):
                writer.writerow([item.get('position', 0), item.get('value', 0)])
        
        return response
    
    return redirect('index')