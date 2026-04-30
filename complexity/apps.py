from django.apps import AppConfig


class ComplexityConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'complexity'
    verbose_name = 'Анализ сложности генома'