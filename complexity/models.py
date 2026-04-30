from django.db import models


class AnalysisResult(models.Model):
    """Модель для сохранения результатов анализа (опционально)"""
    
    METHOD_CHOICES = [
        ('linguistic', 'Лингвистическая сложность'),
        ('shannon', 'Энтропия Шеннона (1-й порядок)'),
        ('shannon2', 'Энтропия Шеннона (2-й порядок)'),
        ('monomer', 'Мера мономеров'),
        ('lz', 'LZ-сложность'),
        ('gc', 'GC-содержание'),
    ]
    
    created_at = models.DateTimeField(auto_now_add=True, verbose_name='Дата создания')
    file_name = models.CharField(max_length=255, verbose_name='Имя файла')
    sequence_length = models.IntegerField(verbose_name='Длина последовательности')
    method = models.CharField(max_length=20, choices=METHOD_CHOICES, verbose_name='Метод')
    window_size = models.IntegerField(verbose_name='Длина окна')
    step = models.IntegerField(verbose_name='Шаг')
    avg_value = models.FloatField(verbose_name='Среднее значение')
    min_value = models.FloatField(verbose_name='Минимум')
    max_value = models.FloatField(verbose_name='Максимум')
    low_regions_count = models.IntegerField(default=0, verbose_name='Количество участков низкой сложности')
    
    class Meta:
        verbose_name = 'Результат анализа'
        verbose_name_plural = 'Результаты анализа'
        ordering = ['-created_at']
    
    def __str__(self):
        return f"{self.file_name} - {self.get_method_display()} - {self.created_at.strftime('%Y-%m-%d %H:%M')}"