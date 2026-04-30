import re
import math

# ==================== БАЗОВЫЕ ФУНКЦИИ ====================

VALID_NUCLEOTIDES = set('ATGC')
AMBIGUOUS_SYMBOLS = set('NRYSWKMBDHV')


def validate_sequence(seq):
    """Проверка, что последовательность состоит только из A, T, G, C"""
    seq_upper = seq.upper()
    
    if not re.fullmatch(r'[ATGC]+', seq_upper):
        invalid_chars = set(seq_upper) - VALID_NUCLEOTIDES
        ambiguous_chars = [c for c in invalid_chars if c in AMBIGUOUS_SYMBOLS]
        
        message = ""
        if ambiguous_chars:
            message += f"Обнаружены консенсусные/неопределённые символы: {set(ambiguous_chars)}. "
            message += "Пожалуйста, используйте последовательность только из A, T, G, C."
        
        return False, list(set(ambiguous_chars)), message
    
    return True, [], ""


# ==================== 1. ЛИНГВИСТИЧЕСКАЯ СЛОЖНОСТЬ ====================

def count_ling_complexity(window_seq, alphabet_size=4):
    """Расчет лингвистической сложности для одного окна"""
    n = len(window_seq)
    if n == 0:
        return 0.0
    
    denominator = 0
    numerator = 0
    
    for i in range(1, n + 1):
        max_possible = min(alphabet_size ** i, n - i + 1)
        denominator += max_possible
        
        unique_words = set()
        for j in range(n - i + 1):
            unique_words.add(window_seq[j:j+i])
        numerator += len(unique_words)
    
    return numerator / denominator if denominator > 0 else 0.0


def calculate_linguistic_complexity(sequence, window_size=50, step=100, alphabet_size=4):
    """Расчет профиля лингвистической сложности"""
    seq_upper = sequence.upper()
    n = len(seq_upper)
    results = []
    
    for start in range(0, n - window_size + 1, step):
        window = seq_upper[start:start + window_size]
        complexity = count_ling_complexity(window, alphabet_size)
        center_pos = start + window_size // 2
        results.append((center_pos, complexity))
    
    if (n - window_size) % step != 0 and n >= window_size:
        window = seq_upper[n - window_size:n]
        complexity = count_ling_complexity(window, alphabet_size)
        center_pos = n - window_size // 2
        results.append((center_pos, complexity))
    
    return results


# ==================== 2. ЭНТРОПИЯ ШЕННОНА ====================

def shannon_entropy(window_seq, alphabet_size=4):
    """Расчет энтропии Шеннона"""
    if not window_seq:
        return 0.0
    
    n = len(window_seq)
    freq = {}
    for char in window_seq.upper():
        freq[char] = freq.get(char, 0) + 1
    
    entropy = 0.0
    for count in freq.values():
        p = count / n
        if p > 0:
            entropy -= p * math.log(p) / math.log(alphabet_size)
    
    return entropy


def shannon_entropy_order2(window_seq, alphabet_size=4):
    """Энтропия второго порядка (учитывает динуклеотиды)"""
    if len(window_seq) < 2:
        return shannon_entropy(window_seq, alphabet_size)
    
    dinuc_freq = {}
    for i in range(len(window_seq) - 1):
        dinuc = window_seq[i:i+2].upper()
        dinuc_freq[dinuc] = dinuc_freq.get(dinuc, 0) + 1
    
    entropy = 0.0
    total = len(window_seq) - 1
    for count in dinuc_freq.values():
        p = count / total
        if p > 0:
            entropy -= p * math.log(p) / math.log(alphabet_size ** 2)
    
    return entropy


def calculate_entropy_profile(sequence, window_size=50, step=100, entropy_order=1, alphabet_size=4):
    """Расчет профиля энтропии Шеннона"""
    seq_upper = sequence.upper()
    n = len(seq_upper)
    results = []
    
    for start in range(0, n - window_size + 1, step):
        window = seq_upper[start:start + window_size]
        if entropy_order == 1:
            entropy = shannon_entropy(window, alphabet_size)
        else:
            entropy = shannon_entropy_order2(window, alphabet_size)
        center_pos = start + window_size // 2
        results.append((center_pos, entropy))
    
    if (n - window_size) % step != 0 and n >= window_size:
        window = seq_upper[n - window_size:n]
        if entropy_order == 1:
            entropy = shannon_entropy(window, alphabet_size)
        else:
            entropy = shannon_entropy_order2(window, alphabet_size)
        center_pos = n - window_size // 2
        results.append((center_pos, entropy))
    
    return results


# ==================== 3. МЕРА МОНОМЕРОВ ====================

def monomer_measure(window_seq):
    """Мера равномерности нуклеотидного состава"""
    if len(window_seq) < 2:
        return 0.0
    
    seq = window_seq.upper()
    changes = 0
    for i in range(len(seq) - 1):
        if seq[i] != seq[i + 1]:
            changes += 1
    
    max_changes = len(seq) - 1
    return changes / max_changes if max_changes > 0 else 0.0


def calculate_monomer_profile(sequence, window_size=50, step=100):
    """Расчет профиля меры мономеров"""
    seq_upper = sequence.upper()
    n = len(seq_upper)
    results = []
    
    for start in range(0, n - window_size + 1, step):
        window = seq_upper[start:start + window_size]
        measure = monomer_measure(window)
        center_pos = start + window_size // 2
        results.append((center_pos, measure))
    
    if (n - window_size) % step != 0 and n >= window_size:
        window = seq_upper[n - window_size:n]
        measure = monomer_measure(window)
        center_pos = n - window_size // 2
        results.append((center_pos, measure))
    
    return results


# ==================== 4. LZ-СЛОЖНОСТЬ (ИСПРАВЛЕННАЯ) ====================

COMPLEMENT_MAP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def lz_complexity_ratio(window_seq, consider_reverse=False, consider_complement=False):
    """
    Нормированная LZ-сложность (алгоритм Лемпеля-Зива)
    Возвращает значение от 0 до 1
    """
    seq = window_seq.upper()
    n = len(seq)
    if n == 0:
        return 0.0
    
    dictionary = set()
    i = 0
    fragments = 0
    
    while i < n:
        # Ограничиваем поиск для скорости (не более 100 символов)
        max_len = min(100, n - i)
        found = False
        
        for length in range(max_len, 0, -1):
            fragment = seq[i:i+length]
            
            # Проверка прямого повтора
            if fragment in dictionary:
                dictionary.add(fragment)
                i += length
                fragments += 1
                found = True
                break
            
            # Проверка инвертированного повтора
            if consider_reverse:
                reverse_fragment = fragment[::-1]
                if reverse_fragment in dictionary:
                    dictionary.add(fragment)
                    i += length
                    fragments += 1
                    found = True
                    break
            
            # Проверка комплементарного повтора
            if consider_complement:
                complement_fragment = ''.join(COMPLEMENT_MAP.get(c, c) for c in fragment)
                if complement_fragment in dictionary:
                    dictionary.add(fragment)
                    i += length
                    fragments += 1
                    found = True
                    break
                
                # Комплементарно-инвертированный
                rev_comp_fragment = ''.join(COMPLEMENT_MAP.get(c, c) for c in fragment[::-1])
                if rev_comp_fragment in dictionary:
                    dictionary.add(fragment)
                    i += length
                    fragments += 1
                    found = True
                    break
        
        # Если не нашли совпадений, добавляем текущий символ
        if not found:
            dictionary.add(seq[i])
            i += 1
            fragments += 1
    
    # Нормировка: максимальное количество фрагментов = длина последовательности
    return fragments / n


def calculate_lz_profile(sequence, window_size=50, step=100, consider_reverse=False, consider_complement=False):
    """Расчет профиля LZ-сложности"""
    seq_upper = sequence.upper()
    n = len(seq_upper)
    results = []
    
    print(f"Расчёт LZ-профиля: длина генома {n}, окно {window_size}, шаг {step}")
    
    for start in range(0, n - window_size + 1, step):
        window = seq_upper[start:start + window_size]
        lz_value = lz_complexity_ratio(window, consider_reverse, consider_complement)
        center_pos = start + window_size // 2
        results.append((center_pos, lz_value))
    
    if (n - window_size) % step != 0 and n >= window_size:
        window = seq_upper[n - window_size:n]
        lz_value = lz_complexity_ratio(window, consider_reverse, consider_complement)
        center_pos = n - window_size // 2
        results.append((center_pos, lz_value))
    
    print(f"LZ-профиль рассчитан: {len(results)} точек")
    return results


# ==================== 5. GC-СОДЕРЖАНИЕ ====================

def gc_content(window_seq):
    """GC-содержание"""
    if len(window_seq) == 0:
        return 0.0
    
    seq = window_seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    return gc_count / len(seq)


def calculate_gc_profile(sequence, window_size=50, step=100):
    """Расчет профиля GC-содержания"""
    seq_upper = sequence.upper()
    n = len(seq_upper)
    results = []
    
    for start in range(0, n - window_size + 1, step):
        window = seq_upper[start:start + window_size]
        gc = gc_content(window)
        center_pos = start + window_size // 2
        results.append((center_pos, gc))
    
    if (n - window_size) % step != 0 and n >= window_size:
        window = seq_upper[n - window_size:n]
        gc = gc_content(window)
        center_pos = n - window_size // 2
        results.append((center_pos, gc))
    
    return results


# ==================== ПОИСК ТАНДЕМНЫХ ПОВТОРОВ ====================

def find_short_tandem_repeats(sequence, min_len=2):
    """
    Поиск коротких тандемных повторов (STR) и гомополимеров
    """
    seq_upper = sequence.upper()
    repeats = []
    
    i = 0
    while i < len(seq_upper):
        # Гомополимеры
        char = seq_upper[i]
        run_length = 1
        while i + run_length < len(seq_upper) and seq_upper[i + run_length] == char:
            run_length += 1
        
        if run_length >= min_len:
            repeats.append({
                'start': i,
                'end': i + run_length,
                'motif': char,
                'length': 1,
                'repeat_count': run_length,
                'type': 'homopolymer'
            })
            i += run_length
            continue
        
        # Динуклеотидные повторы
        if i + 2 < len(seq_upper):
            dinuc = seq_upper[i:i+2]
            if dinuc[0] != dinuc[1]:
                count = 1
                pos = i + 2
                while pos + 2 <= len(seq_upper) and seq_upper[pos:pos+2] == dinuc:
                    count += 1
                    pos += 2
                if count >= 2:
                    repeats.append({
                        'start': i,
                        'end': pos,
                        'motif': dinuc,
                        'length': 2,
                        'repeat_count': count,
                        'type': 'dinucleotide'
                    })
                    i = pos
                    continue
        
        # Тринуклеотидные повторы
        if i + 3 < len(seq_upper):
            trinuc = seq_upper[i:i+3]
            if len(set(trinuc)) > 1:
                count = 1
                pos = i + 3
                while pos + 3 <= len(seq_upper) and seq_upper[pos:pos+3] == trinuc:
                    count += 1
                    pos += 3
                if count >= 2:
                    repeats.append({
                        'start': i,
                        'end': pos,
                        'motif': trinuc,
                        'length': 3,
                        'repeat_count': count,
                        'type': 'trinucleotide'
                    })
                    i = pos
                    continue
        
        i += 1
    
    # Сортируем по размеру
    repeats.sort(key=lambda x: (x['end'] - x['start']), reverse=True)
    
    return repeats


# ==================== ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ====================

def parse_bed_file(bed_file):
    """Парсинг BED файла"""
    genes = []
    
    if hasattr(bed_file, 'read'):
        content = bed_file.read()
        if isinstance(content, bytes):
            content = content.decode('utf-8')
    elif isinstance(bed_file, bytes):
        content = bed_file.decode('utf-8')
    elif isinstance(bed_file, str):
        content = bed_file
    else:
        return genes
    
    lines = content.split('\n')
    
    for line in lines:
        line = line.strip()
        if line and not line.startswith('#'):
            parts = line.split('\t')
            if len(parts) >= 4:
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                    name = parts[3]
                    strand = parts[5] if len(parts) > 5 else '+'
                    genes.append({
                        'start': start,
                        'end': end,
                        'name': name,
                        'strand': strand,
                        'length': end - start
                    })
                except ValueError:
                    continue
    
    return genes


def get_statistics(profile):
    """Получение статистики по профилю"""
    values = [v for _, v in profile]
    if not values:
        return {'avg': 0, 'min': 0, 'max': 0}
    return {
        'avg': sum(values) / len(values),
        'min': min(values),
        'max': max(values)
    }


def find_local_minima(profile, window=5, threshold=0.7):
    """Поиск локальных минимумов сложности"""
    minima = []
    
    for i in range(len(profile)):
        pos, val = profile[i]
        
        if val >= threshold:
            continue
        
        is_minimum = True
        left_start = max(0, i - window)
        right_end = min(len(profile), i + window + 1)
        
        for j in range(left_start, right_end):
            if j != i and profile[j][1] < val:
                is_minimum = False
                break
        
        if is_minimum:
            minima.append((pos, val, i))
    
    return minima


def get_sequence_region(sequence, center_pos, window=30):
    """Получить фрагмент последовательности вокруг позиции"""
    start = max(0, center_pos - window)
    end = min(len(sequence), center_pos + window)
    region = sequence[start:end]
    
    center_in_region = center_pos - start
    if 0 <= center_in_region < len(region):
        highlighted = region[:center_in_region] + f"[{region[center_in_region]}]" + region[center_in_region+1:]
    else:
        highlighted = region
    
    return highlighted, start, end