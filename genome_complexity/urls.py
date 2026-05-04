from django.urls import path
from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('about/', views.about, name='about'),
    path('analyze/', views.analyze, name='analyze'),
    path('export/', views.export_results, name='export_results'),
]