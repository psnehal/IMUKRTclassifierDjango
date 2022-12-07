from django.db import models

# Create your models here.


class Jobinput(models.Model):
    uuid = models.TextField()
    pcaparam=models.TextField()
    logparam=models.TextField()
    infilepath =models.TextField()


