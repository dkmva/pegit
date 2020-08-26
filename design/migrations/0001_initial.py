# Generated by Django 3.1 on 2020-08-26 06:33

import design.models
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Gene',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('gene_id', models.CharField(db_index=True, max_length=200)),
                ('name', models.CharField(db_index=True, max_length=50)),
                ('chromosome', models.CharField(max_length=200)),
                ('strand', models.CharField(max_length=1)),
                ('start', models.PositiveIntegerField()),
                ('end', models.PositiveIntegerField()),
                ('gene_type', models.CharField(max_length=200)),
                ('source', models.CharField(max_length=200)),
            ],
            options={
                'ordering': ('id',),
            },
            bases=(models.Model, design.models.SequenceObjectMixin),
        ),
        migrations.CreateModel(
            name='Organism',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=200)),
                ('assembly', models.CharField(max_length=200)),
                ('annotation', models.CharField(max_length=200)),
                ('source', models.CharField(max_length=200)),
                ('codon_table', models.CharField(max_length=200, null=True)),
                ('sequence_search', models.CharField(max_length=200)),
                ('scaffolds', models.TextField(default='')),
            ],
            options={
                'ordering': ('id',),
            },
        ),
        migrations.CreateModel(
            name='Transcript',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('transcript_id', models.CharField(db_index=True, max_length=200)),
                ('name', models.CharField(db_index=True, max_length=200)),
                ('start', models.PositiveIntegerField()),
                ('end', models.PositiveIntegerField()),
                ('transcript_type', models.CharField(max_length=200)),
                ('source', models.CharField(max_length=200)),
                ('gene', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='transcripts', to='design.gene')),
            ],
            options={
                'ordering': ('transcript_id',),
            },
            bases=(models.Model, design.models.SequenceObjectMixin),
        ),
        migrations.AddField(
            model_name='gene',
            name='organism',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='genes', to='design.organism'),
        ),
        migrations.CreateModel(
            name='Exon',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('start', models.PositiveIntegerField()),
                ('end', models.PositiveIntegerField()),
                ('transcript', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='exons', to='design.transcript')),
            ],
        ),
        migrations.CreateModel(
            name='CodingSequence',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('start', models.PositiveIntegerField()),
                ('end', models.PositiveIntegerField()),
                ('transcript', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='coding_sequences', to='design.transcript')),
            ],
        ),
    ]