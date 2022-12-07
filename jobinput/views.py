import mimetypes

from django.http import HttpResponse
from django.http import JsonResponse

from django.shortcuts import render
from .forms import PostForm
from django.core.files.storage import FileSystemStorage
import os
from django.conf import settings
import uuid
import subprocess
from django.core.files import File



# Create your views here.

#IMU_KRT_classifier.py [-h] -dir DIR -PCA PCA -log2cpmmatrix INPUT

def home(request,*args,**kwargs):
    print("its in the home")
    print(args, kwargs)
    form = PostForm(request.POST)
    #return HttpResponse("<h1>Hello World</h1>") # string of HTML code
    return render (request, "home.html",{"form": form })


def heatmap(request,*args,**kwargs):

    return render (request, "heatmap.html",{})


def download_file(request,imageid):

    if imageid != '':
        print("filepath is "+imageid)

        #https://stackoverflow.com/questions/63228094/how-can-i-create-download-link-using-django
        filepath = '/Users/snehalpatil/Documents/GithubProjects/ShitingProject/tumourclass/cpm_input_classifier.csv'

        # Open the file for reading content
        path = open('/Users/snehalpatil/Documents/GithubProjects/ShitingProject/tumourclass/cpm_input_classifier.csv', "r")
        # Set the mime type
        mime_type, _ = mimetypes.guess_type(filepath)
        print (mime_type)
        # Set the return value of the HttpResponse
        response = HttpResponse(path, content_type=mime_type)
        # Set the HTTP header for sending to browser
        response['Content-Disposition'] = "attachment; filename=%s" % filepath
        return response
    else:
        # Load the template
        return render(request, 'index.html')




form = PostForm()

def new(request):

    return render(request, "posts/new.html", )

FILE_UPLOAD_DIR = '/home/imran/uploads'



# ajax_posting function
def ajax_posting(request):
    if request.is_ajax():
        first_name = request.POST.get('first_name', None) # getting data from first_name input
        last_name = request.POST.get('last_name', None)  # getting data from last_name input
        if first_name and last_name: #cheking if first_name and last_name have value
            response = {
                'msg':'Your form has been submitted successfully' # response message
            }
            return JsonResponse(response) # return response as JSON




def jobsubmit(request):
#    <QueryDict: {'batcheff': ['no'], 'filet': ['0'], 'logcpm': ['1'], 'csrfmiddlewaretoken': ['MIKpWc5gL8uhQfxTWuCqbYEUGmj7P2aIkbJ4OYpC70c6VTlHrVFjCOvCe6wNJT9u']}>
    print("hi from job submit")
    print(request.POST)



    #ileinput[absolute path],filesep[comma,tab],cpmornot[0,1],
    # batchremovalornot(has to provide raw count)[combat(include FF and cannot be changed),ruvg,not],includeFF18ornot[0,1] \




    logcpm = '1'
    method = 1


    lines = []
    d= {}
    myfile='/Users/snehalpatil/Documents/GithubProjects/ShitingProject/IMUKRTclassifier/demo_data/all_samples_2536.count'
    if request.method == 'POST' and request.FILES['uploadfile']:
    #if logcpm == '1':
        myfile = request.FILES['uploadfile']
        print("myfile")
        print(myfile)

        folder =  str(uuid.uuid4())




        try:
            os.mkdir(os.path.join(settings.MEDIA_ROOT, folder))
        except:
            pass

        file_save_path = os.path.join(settings.MEDIA_ROOT, folder)
        print(file_save_path)
        #fs = FileSystemStorage()
        fs = FileSystemStorage(location=file_save_path , base_url = file_save_path) #defaults to   MEDIA_ROOT
        print(myfile)
        filename = fs.save('all_samples_2536.count', myfile)

        uploaded_file_url = fs.url(filename)
        print("print(uploaded_file_url)")
        print(uploaded_file_url)

    # with open('/Users/snehalpatil/Documents/GithubProjects/ShitingProject/things_give_snehal/output3.txt') as f:
    #     lines = f.readlines()
    # count = 0
    # for line in lines:
    #     count += 1
    #     #print(f'line {count}: {line.split(", ")[0]}')
    #     (key, val) = line.split()
    #     d[key] = val

    #/Users/snehalpatil/Documents/GithubProjects/ShitingProject/tumourclass/jobinput/static/jobinput/R_script/preprocessing.R
# /Users/snehalpatil/Documents/GithubProjects/ShitingProject/tumourclass/media/87e39b5c-2753-4e92-ae7e-e87943e96099/all_samples_2536.count
# 0  not NULL 0  /Users/snehalpatil/Documents/GithubProjects/ShitingProject/things_give_snehal/


        stmt  = "Rscript /Users/snehalpatil/Documents/GithubProjects/ShitingProject/tumourclass/jobinput/static/jobinput/R_script/preprocessing.R "
        stmt+= uploaded_file_url
        stmt+=   " 0  "
        stmt+= "not NULL 0 "
        stmt+= file_save_path

        print("stmt :"+stmt)

        res = subprocess.call(stmt,shell = True)
        print("*************************************")
        print(res)
        django_file =''

        outputfile = file_save_path + '/cpm_input_classifier.csv'
        if(res == 0):
            some_file = open(outputfile, "r")
            django_file = File(some_file)


    filepath = '/Users/snehalpatil/Documents/GithubProjects/ShitingProject/tumourclass/cpm_input_classifier.csv'
    print("django_file",django_file)
    # Open the file for reading content
    path = open('/Users/snehalpatil/Documents/GithubProjects/ShitingProject/tumourclass/cpm_input_classifier.csv', "r")
    # Set the mime type
    mime_type = mimetypes.guess_type(filepath)
    # Set the return value of the HttpResponse

    #response = HttpResponse(path, content_type=mime_type)
    # Set the HTTP header for sending to browser
    #response['Content-Disposition'] = "attachment; filename=%s" % myfile
    # Return the response value
    #return response
    #return render (request, "jobsubmit.html",{"fname":django_file })# string of HTML code
    #return outputfile

    return JsonResponse({"outputfile": outputfile}, status=200)



    #return render (request, "jobsubmit.html",{"list": lines})# string of HTML code


def clustergramheatmap(request):
    print ("testng clustergram")
    return render (request, "clustergramheatmap.html",{"form": "tes" })




def runclassifier(request):

    print(request.POST)
    lines = []
    d= {}
    import glob
    import pandas as pd


    cols = [0, 1, 2] # add more columns here

    df = pd.DataFrame()
    arr = pd.read_csv("/Users/snehalpatil/Documents/GithubProjects/ShitingProject/IMUKRTclassifier/demo_result/predict_result.txt", sep='\t', header=None, usecols=cols)
    for ind in arr.index:
      print(arr[0][ind],arr[1][ind], arr[2][ind])
    # with open('/Users/snehalpatil/Documents/GithubProjects/ShitingProject/IMUKRTclassifier/demo_result/predict_result.txt') as f:
    #     lines = f.readlines()
    #     count = 0
    #     for line in lines:
    #         count += 1
    #         print(f'line {count}: {line.split()}')
    #         (key, val) = line.split()
    #         d[key] = val
    return render(request, "runclassifier.html", {"arr": arr} )


