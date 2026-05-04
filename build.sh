#!/bin/bash

# Устанавливаем зависимости
pip install -r requirements.txt

# Собираем статические файлы
python manage.py collectstatic --noinput

# Выполняем миграции
python manage.py migrate