from django.conf.urls import url, include
from django.urls import path, re_path
from django.views.generic import RedirectView
from rest_framework import routers

from . import views

router = routers.DefaultRouter()
router.register(r'organisms', views.OrganismViewSet, basename='organism')
router.register(r'genes', views.GeneViewSet, basename='gene')
router.register(r'transcripts', views.TranscriptViewSet, basename='transcript')
router.register(r'edits', views.EditViewSet, basename='edit')
router.register(r'jobs', views.JobViewSet, basename='job')
router.register(r'clinvar', views.ClinvarSearch, basename='clinvar')


favicon_view = RedirectView.as_view(url='/static/favicon.ico', permanent=True)
favicon32_view = RedirectView.as_view(url='/static/favicon-32x32.png', permanent=True)
favicon16_view = RedirectView.as_view(url='/static/favicon-16x16.png', permanent=True)
apple_touch_view = RedirectView.as_view(url='/static/apple-touch-icon.png', permanent=True)
gencode_example_file = RedirectView.as_view(url='/static/example_edits_gencode_release_33', permanent=True)
clinvar_example_file = RedirectView.as_view(url='/static/example_edits_clinvar', permanent=True)


urlpatterns = [
    path('', views.index, name='index'),
    path('instructions', views.index, name='index'),
    path('about', views.index, name='index'),
    path('contact', views.index, name='index'),
    path('faq', views.index, name='index'),
    path('customsequence/', views.index, name='index'),
    re_path('job/.*', views.index, name='index'),
    re_path('transcript/.*', views.index, name='index'),
    re_path('gene/.*', views.index, name='index'),
    url(r'^favicon\.ico$', favicon_view),
    url(r'^favicon-32x32\.png$', favicon32_view),
    url(r'^favicon-16x16\.png$', favicon16_view),
    url(r'^apple-touch-icon\.png$', apple_touch_view),
    url(r'^example_edits_gencode_release_33$', gencode_example_file),
    url(r'^example_edits_clinvar$', clinvar_example_file),
    path(r'api/', include(router.urls)),
    re_path('static/.*', include('django.contrib.staticfiles.urls')),
    #path('api/clinvar/', views.ClinvarSearch.as_view()),
    #re_path(r'[^api/]', views.index, name='index'),
]

