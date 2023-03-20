import mimetypes

from django.http import HttpResponse
from django.http import JsonResponse
from django.shortcuts import redirect

from django.shortcuts import render


from .forms import PostForm
from django.core.files.storage import FileSystemStorage
import os
from django.conf import settings
import uuid
import subprocess
from django.core.files import File

from .python_script.IMU_KRT_classifier import runfile

import glob
import pandas as pd




# ...
#
# try:
#     SECRET_KEY = os.environ["SECRET_KEY"]
# except KeyError as e:
#     raise RuntimeError("Could not find a SECRET_KEY in environment") from e



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


def about(request,*args,**kwargs):

    return render (request, "about.html",{})


def download_file(uuidno):



    if uuidno != '':
        print("filepath is "+uuidno)

        #https://stackoverflow.com/questions/63228094/how-can-i-create-download-link-using-django
        filepath = '/Users/snehalpatil/Documents/GithubProjects/ShitingProject/tumourclass/media/'+uuid+'/cpm_input_classifier.csv'

        # Open the file for reading content
        path = open(filepath, "r")
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
        return render('index.html')




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
    print("hi from job submit",request.POST)
    print(request.POST.get('batcheff'))

    lines = []
    d= {}
    res=0
    #myfile='/Users/snehalpatil/Documents/GithubProjects/ShitingProject/IMUKRTclassifier/demo_data/all_samples_2536.count'
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
        filename = fs.save(myfile.name, myfile)
        uploaded_file_url = fs.url(filename)
        print(uploaded_file_url)




# message('usage: fileinput[absolute path],cpmornot[0,1],batchremovalornot(has to provide raw count)[combat,ruvg,not],batch_effect_file[either input a NULL or a directory to files store batcheffect],includeFF18ornot[0,1],output_directory[absolute path],folder_store_required_information[abosulte path]')
#fileinput,=uploaded file
# cpmornot =Indicate file type
# batchremovalornot(has to provide raw count)[combat,ruvg,not],=Do you have multiple batches? logic batchremovalornot
# batch_effect_file[either input a NULL or a directory to files store batcheffect],=filepath2
# includeFF18ornot[0,1]==Include pre-loaded additional samples? preloads
# output_directory[absolute path],
# folder_store_required_information[abosulte path]')


        batcheff =request.POST.get('batcheff')
        cpmornot=request.POST.get('logcpm')
        includeFF18ornot=request.POST.get('preloads')
        batchfilestatus=request.POST.get('batchfilestatus')

        print("settings.path",settings.RSCRIPT_PATH)
        rscriptfilepath = os.path.join(settings.RSCRIPT_PATH, "R_script/")+"preprocessing.R"
        print("rscriptfilepath",rscriptfilepath)

  #batchremovalornot logic
        if (includeFF18ornot == 1 and batcheff == 1 ):
            print("yes")
            batchremovalornot = "combat"
        elif(includeFF18ornot == 0 and batcheff == 1 ):
            batchremovalornot = "combat"
        else:
            batchremovalornot = "not"


        if (batchfilestatus == "yes"):
            print("got the batch file")
            batch_effect_file = request.FILES['batchfile']
            print("batch_effect_file")
            print(batch_effect_file)
        else:
            print("batch file is Null")
            batch_effect_file = 'NULL'


        stmt  = "Rscript  "
        stmt+= rscriptfilepath + " "
        stmt+= uploaded_file_url+ " "
        stmt+= cpmornot+ " "
        stmt+= batchremovalornot+ " "
        stmt+= batch_effect_file+ " "
        stmt+= includeFF18ornot+ " "
        stmt+= file_save_path+ " "
        stmt+= os.path.join(settings.RSCRIPT_PATH, "R_script/")
        print("stmt :"+stmt)

        res = subprocess.call(stmt,shell = True)
        print("*************************************")
        print(res)
        django_file =''

        outputfile = file_save_path + '/cpm_input_classifier.csv'
        if(res == 0):
            some_file = open(outputfile, "r")
            django_file = File(some_file)
            #print("django_file",django_file)

        else:
            print("its in th else loop")
            response = redirect('DisplayError')
    print(res)
    return JsonResponse({"foldername": folder,"res":res})


def DisplayError(request):
    print("reached to the error")
    return render (request, "DisplayError.html")


    #return render (request, "jobsubmit.html",{"list": lines})# string of HTML code


def clustergramheatmap(request):
    print ("testng clustergram")
    return render (request, "clustergramheatmap.html",{"form": "tes" })




def runclassifier(request):

    print(request.POST.get('outputfile'))
    lines = []
    d= {}

    folderpath = request.POST.get('outputfile')

    #result file from the
    file_save_path = os.path.join(settings.MEDIA_ROOT, folderpath,"cpm_input_classifier.csv" )

    print(file_save_path)

    pyfilepath = os.path.join(settings.RSCRIPT_PATH, "python_script/")
    print("pyfilepath",pyfilepath)



    #python3 IMU_KRT_classifier.py -dir ./ -PCA 1 -log2cpmmatrix ../demo_data/cpm_input_classifier.csv -output_dir ../demo_result/



    cols = [0, 1, 2] # add more columns here

    df = pd.DataFrame()
    arr = pd.read_csv("/Users/snehalpatil/Documents/GithubProjects/ShitingProject/IMUKRTclassifier/demo_result/predict_result.txt", sep='\t', header=None, usecols=cols)
    # for ind in arr.index:
    #   print(arr[0][ind],arr[1][ind], arr[2][ind])
    # with open('/Users/snehalpatil/Documents/GithubProjects/ShitingProject/IMUKRTclassifier/demo_result/predict_result.txt') as f:
    #     lines = f.readlines()
    #     count = 0
    #     for line in lines:
    #         count += 1
    #         print(f'line {count}: {line.split()}')
    #         (key, val) = line.split()
    #         d[key] = val

    #python3 IMU_KRT_classifier.py -dir ./ -PCA 1 -log2cpmmatrix samples_result/18_plus_FFPE.csv -output_dir ./
    #workdir, pca, input,output_dir,models,gene_set

    # python3 IMU_KRT_classifier.py
    # -dir /Users/snehalpatil/Documents/GithubProjects/ShitingProject/IMUKRTclassifier-lab_version/IMUKRT_classifier_lab_version/python_script
    # -PCA 1
    # -log2cpmmatrix  /Users/snehalpatil/Documents/GithubProjects/ShitingProject/IMUKRTclassifier-lab_version/IMUKRT_classifier_lab_version/demo_data/cpm_input_classifier.csv
    # -output_dir /Users/snehalpatil/Documents/GithubProjects/ShitingProject/tumourclass/media/2d14ba8a-b589-4fe6-a28b-4d654c77ba10/
    # -models rf,knn,gnb,svm
    # -genes 168,960

    models ='f,knn,gnb,svm,elasticnet'
    gene_set ='168,960'
    runfile(pyfilepath,"1","/Users/snehalpatil/Documents/GithubProjects/ShitingProject/IMUKRTclassifier-lab_version/IMUKRT_classifier_lab_version/results/cpm_input_classifier.csv","/Users/snehalpatil/Documents/GithubProjects/ShitingProject/IMUKRTclassifier-lab_version/IMUKRT_classifier_lab_version/results/",None,None)
    return render(request, "runclassifier.html", {"arr": arr} )


