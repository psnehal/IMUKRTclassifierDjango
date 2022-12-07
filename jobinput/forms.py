from django.forms import ModelForm, Textarea
from .models import Jobinput

class PostForm(ModelForm):
    class Meta:
        model = Jobinput
        fields = '__all__'