# Generated by Django 3.1 on 2021-03-01 21:05

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('design', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='gene',
            name='name',
            field=models.CharField(db_index=True, max_length=200),
        ),
    ]