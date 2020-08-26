from django.contrib import admin

# Register your models here.
from .models import Organism, Gene, Transcript, CodingSequence


@admin.register(Organism)
class OrganismAdmin(admin.ModelAdmin):
    pass


admin.site.register(Gene)
admin.site.register(Transcript)
admin.site.register(CodingSequence)


#@admin.register(Job)
#class JobAdmin(admin.ModelAdmin):
    #list_display = ('job_id', 'result')
    #raw_id_fields = ('gene', )
