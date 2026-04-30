from django import forms


METHOD_CHOICES = [
    ('linguistic', 'Лингвистическая сложность (Troyanskaya et al., 2002)'),
    ('shannon', 'Энтропия Шеннона (1-й порядок)'),
    ('shannon2', 'Энтропия Шеннона (2-й порядок)'),
    ('monomer', 'Мера мономеров (равномерность состава)'),
    ('gc', 'GC-содержание'),
]


class FastaUploadForm(forms.Form):
    fasta_file = forms.FileField(
        label='FASTA файл с последовательностью',
        help_text='Файл должен содержать одну нуклеотидную последовательность (A, T, G, C)',
        widget=forms.FileInput(attrs={'accept': '.fasta,.fa,.txt'})
    )
    
    method = forms.ChoiceField(
        choices=METHOD_CHOICES,
        label='Метод анализа сложности',
        initial='linguistic',
        widget=forms.Select(attrs={'class': 'method-select'})
    )
    
    window_size = forms.IntegerField(
        label='Длина окна (нт)',
        min_value=10,
        max_value=1000,
        initial=50,
        help_text='Рекомендуемое значение: 50'
    )
    
    step = forms.IntegerField(
        label='Шаг скольжения (нт)',
        min_value=1,
        max_value=500,
        initial=100,
        help_text='Рекомендуемое значение: 100'
    )
    
    alphabet_size = forms.IntegerField(
        label='Размер алфавита (для энтропии)',
        min_value=2,
        max_value=4,
        initial=4,
        required=False,
        help_text='Для ДНК/РНК используйте 4'
    )
    
    bed_file = forms.FileField(
        label='BED файл с аннотацией генов (опционально)',
        required=False,
        help_text='Файл с координатами генов для отображения на графике',
        widget=forms.FileInput(attrs={'accept': '.bed,.txt'})
    )
    
    def clean_fasta_file(self):
        fasta_file = self.cleaned_data['fasta_file']
        if fasta_file.size > 10 * 1024 * 1024:
            raise forms.ValidationError('Размер файла не должен превышать 10 МБ')
        return fasta_file
    
    def clean_window_size(self):
        window_size = self.cleaned_data['window_size']
        if window_size < 10:
            raise forms.ValidationError('Длина окна должна быть не менее 10 нуклеотидов')
        return window_size
    
    def clean_step(self):
        step = self.cleaned_data['step']
        if step < 1:
            raise forms.ValidationError('Шаг должен быть положительным числом')
        return step