"""csx_unchained URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.9/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.urls import include, re_path
from django.contrib import admin
from django.conf import settings
from django.conf.urls.static import static

import consensx.views

urlpatterns = [
    re_path(r'^admin/', admin.site.urls),
    re_path(r'^$', consensx.views.home),
    re_path(r'selection/(?P<my_id>[a-zA-Z0-9]{6}$)', consensx.views.selection),
    re_path(r'^(?P<my_id>[a-zA-Z0-9]{6}$)', consensx.views.db),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
