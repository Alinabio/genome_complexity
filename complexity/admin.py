from django.contrib import admin
from .models import AnalysisResult

@admin.register(AnalysisResult)
class AnalysisResultAdmin(admin.ModelAdmin):
    list_display = ['file_name', 'method', 'sequence_length', 'avg_value', 'created_at']
    list_filter = ['method', 'created_at']
    search_fields = ['file_name']
    readonly_fields = ['created_at']