"""tumourclass URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""

from django.urls import path

from jobinput.views import ajax_posting, jobsubmit
from . import views


urlpatterns = [
    path('', views.home, name='home'),
    path('jobsubmit/', jobsubmit, name='jobsubmit'),
    path('runclassifier',views.runclassifier, name='runclassifier'),
    path('clustergramheatmap',views.clustergramheatmap, name='clustergramheatmap'),
    path('heatmap', views.heatmap, name='heatmap'),
    path('downloadfile/<str:imageid>/', views.download_file, name='download_file'),
    path('ajax-posting/', ajax_posting, name='ajax_posting'),



]
