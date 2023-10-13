import mimetypes
import zipfile
from io import BytesIO

from zipfile import ZipFile

from django.http import HttpResponse
from django.http import FileResponse
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
from clustergrammer import Network

import glob
import pandas as pd
import json

from django.shortcuts import render

from urllib.parse import quote




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

   #print("pcaf", pcaf)
   pcaf = '/Users/snehalpatil/Documents/GithubProjects/ShitingProject/tumourclass/media/be8f318a-604d-47f6-8c2e-5d49c7283e78/mult_view.json'
   with open(pcaf) as jsonFile:
       pcadata = json.load(jsonFile)
       foldername= '96555372-c8ac-4a2b-81c3-ab357929c1dc'

       encoded_foldername = quote(foldername)
       print(pcadata);


   return render (request, "heatmap.html",{"pcadata": pcadata})


def about(request,*args,**kwargs):

    return render (request, "about.html",{})

def download_result(request):
    foldername = request.GET.get('foldername')
    print(foldername)
    folderpath= os.path.join(settings.MEDIA_ROOT, foldername)
    file_blueprint = BytesIO()
    zip_file = ZipFile(file_blueprint, 'a')

#create the full path to the folder you want to download the files from
    filenames = ['cpm_input_classifier_no_batchre.csv','cpm_input_classifier.csv','heatmap_data.txt','predict_result.txt','PCA_classification_result_of_IMU_KRT.png','heatmap_for_clustering_result.png']

    for filename in os.listdir(folderpath):
        try:
            #characterize the from path and the destination path as first and second argument
            if filename in filenames:
                zip_file.write(os.path.join(folderpath + "/" + filename), arcname=filename)
        except Exception as e:
            print(e)
    zip_file.close()
    response = HttpResponse(file_blueprint.getvalue(), content_type = "application/x-zip-compressed")
    response["Content-Disposition"] = "attachment; filename= subtype_classified_results.zip"
    return response

def downloadgeneset(request):
    filename = request.GET.get('filename')
    file_save_path = os.path.join(settings.STATIC_ROOT, filename )
    # Set the mime type
    mime_type, _ = mimetypes.guess_type(file_save_path)
    print (filename)
    # Set the return value of the HttpResponse

    if(filename=="sample.csv"):

        response = HttpResponse(file_save_path, content_type='csv')
        # Set the HTTP header for sending to browser
        response['Content-Disposition'] = "attachment; filename=%s" , filename
    elif(filename=="geneset.xlsx"):
        response = FileResponse(open(file_save_path, 'rb'), content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
        response['Content-Disposition'] = 'attachment; filename="TrainingGeneList.xlsx"'
    else:
        print("inside predict loop loop ")
        file_save_path = os.path.join(settings.MEDIA_ROOT, filename,"predict_result.txt" )
        print(file_save_path)
        response = FileResponse(open(file_save_path, 'rb'), content_type='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
        response['Content-Disposition'] = 'attachment; filename="predict_result.txt"'




    return response




def download_file(request):

    uuidno = request.GET.get('uuid')
    batchrem= request.GET.get('batchrem')
    if uuidno != '':
        print(uuidno)
        if(batchrem =='combat'):
            file_save_path = os.path.join(settings.MEDIA_ROOT, uuidno,"cpm_input_classifier.csv" )

        else:
            file_save_path = os.path.join(settings.MEDIA_ROOT, uuidno,"cpm_input_classifier_no_batchre.csv" )


        #https://stackoverflow.com/questions/63228094/how-can-i-create-download-link-using-django
        #filepath = '/Users/snehalpatil/Documents/GithubProjects/ShitingProject/tumourclass/media/'+str(uuidno)+'/cpm_input_classifier.csv'

        # Open the file for reading content
        path = open(file_save_path, "r")
        filename = os.path.basename(file_save_path)
        # Set the mime type
        mime_type, _ = mimetypes.guess_type(file_save_path)
        print (mime_type)
        # Set the return value of the HttpResponse
        response = HttpResponse(path, content_type=mime_type)
        # Set the HTTP header for sending to browser
        response['Content-Disposition'] = "attachment; filename=%s" , filename
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

        #print("settings.path",settings.RSCRIPT_PATH)
        rscriptfilepath = os.path.join(settings.RSCRIPT_PATH, "R_script/")+"preprocessing1.R"
        print("batcheff",type(batcheff))
        print("includeFF18ornot",includeFF18ornot)

  #batchremovalornot logic
        if (int(batcheff) == 1) or int(includeFF18ornot) == 1:
            print("yes")
            batchremovalornot = "combat"
            #if they select batcheffect == 0 and includeFF18ornot == 0 then batchremovalornot= not
        else:
            batchremovalornot = "not"


        if (batchfilestatus == "yes"):
            print("got the batch file")
            batch_effect_file = request.FILES['batchfile']
            print("batch_effect_file")
            print(batch_effect_file)
            if request.method == 'POST' and request.FILES['batchfile']:
                #if logcpm == '1':
                batch_effect_file = request.FILES['batchfile']
                print("batch_effect_file")
                print(batch_effect_file)
                batchfile = fs.save(batch_effect_file.name, batch_effect_file)
                uploaded_file_urlb = fs.url(batchfile)
                print("uploaded_file_urlb", uploaded_file_urlb)
            else:
                print("batch file is Null")
                uploaded_file_urlb= 'NULL'
                #batch_effect_file = 'NULL'
        else:
            uploaded_file_urlb= 'NULL'


        stmt  = "Rscript  "
        stmt+= rscriptfilepath + " "
        stmt+= uploaded_file_url+ " "
        stmt+= cpmornot+ " "
        stmt+= batchremovalornot+ " "
        stmt+= uploaded_file_urlb+ " "
        stmt+= includeFF18ornot+ " "
        stmt+= file_save_path+ " "
        stmt+= os.path.join(settings.RSCRIPT_PATH, "R_script/")
        print("stmt :"+stmt)



        res = subprocess.call(stmt,shell = True)
        print("*************************************")
        print(res)
        django_file =''



        pcaf = file_save_path + '/PCA_related.json'
        outputfile =''
        outputfilename=''
        pcadata=''
        finalstringno=''
        if(res ==0):
            print("res is 0")
            if(res == 0 and batchremovalornot == 'not' ):
                outputfile = file_save_path + '/cpm_input_classifier_no_batchre.csv'
                print("its in not batchremovalornot loop ")
                some_file = open(outputfile, "r")
                django_file = File(some_file)
                #print("pcaf", pcaf)
                pcafile = File(pcaf)
                with open(pcaf) as jsonFile:
                    pcadata = json.load(jsonFile)
            elif(res == 0 and batchremovalornot == 'combat' ):
                outputfile =file_save_path  + '/raw_count_input.csv'
                print("its in batchremovalornot yes loop")
                some_file = open(outputfile, "r")
                django_file = File(some_file)
                #print("pcaf", pcaf)
                pcafile = File(pcaf)
                with open(pcaf) as jsonFile:
                    pcadata = json.load(jsonFile)
            keytotal= 0
            valuetotal=0
            for (k, v) in pcadata.items():
                print("Key: " + k)
                if(k != "PCA" and k !="warning"):
                    keytotal = keytotal+int(k)
                    valuetotal= valuetotal+v
                    print("Value: " + str(v))

            print("keytotal",keytotal)
            print("valuetotal",valuetotal)
            finalno=  round((valuetotal*100)/keytotal,2)
            print("finalno: " + str(finalno))
            finalstringno = str(valuetotal)+" out of "+str(keytotal)

        else:
            print("res is 1")
            print("its in th else loop")
            response = redirect('DisplayError')






    return JsonResponse({"foldername": folder,"res":res,"pcajson":pcadata,"batchremovalornot":batchremovalornot,outputfile:outputfile,"finalperc":finalstringno})


def runstepone(request):
    #    <QueryDict: {'batcheff': ['no'], 'filet': ['0'], 'logcpm': ['1'], 'csrfmiddlewaretoken': ['MIKpWc5gL8uhQfxTWuCqbYEUGmj7P2aIkbJ4OYpC70c6VTlHrVFjCOvCe6wNJT9u']}>
    print("hi from runstepone",request.POST)

    lines = []
    d= {}
    res=0
    foldername = request.POST.get('foldername')
    file_save_path = os.path.join(settings.MEDIA_ROOT, foldername)
    filename =  file_save_path+'/raw_count_input.csv'
    batcheff =request.POST.get('batcheff')
    cpmornot=request.POST.get('logcpm')
    includeFF18ornot=request.POST.get('preloads')
    batchfilestatus=request.POST.get('batchfilestatus')
    rscriptfilepath = os.path.join(settings.RSCRIPT_PATH, "R_script/")+"preprocessing2.R"
    batch_file_name =request.POST.get('batchfile')
    batchremovalornot=request.POST.get('batchremovalornot')
    if (batchfilestatus == "yes"):

        batchfile = file_save_path+"/"+batch_file_name
    else:

        batchfile= 'NULL'
            #batch_effect_file = 'NULL'



    stmt  = "Rscript  "
    stmt+= rscriptfilepath + " "
    stmt+= filename+ " "
    stmt+= cpmornot+ " "
    stmt+= batchremovalornot+ " "
    stmt+= batchfile+ " "
    stmt+= includeFF18ornot+ " "
    stmt+= file_save_path+ " "
    stmt+= os.path.join(settings.RSCRIPT_PATH, "R_script/")
    print("stmt :"+stmt)



    res = subprocess.call(stmt,shell = True)
    print("*************************************")
    print(res)
    django_file =''



    pcaf = file_save_path + '/PCA_related.json'
    outputfile =''
    if(res == 0 and batchremovalornot == 'not' ):
        outputfile = file_save_path + '/cpm_input_classifier_batch.csv'
        print("its in not batchremovalornot loop ")
        some_file = open(outputfile, "r")
        django_file = File(some_file)
        #print("pcaf", pcaf)
        pcafile = File(pcaf)
        with open(pcaf) as jsonFile:
            pcadata = json.load(jsonFile)
    elif(res == 0 and batchremovalornot == 'combat' ):
        outputfile = file_save_path + '/raw_count_input.csv'
        print("its in batchremovalornot yes loop")
        some_file = open(outputfile, "r")
        django_file = File(some_file)
        #print("pcaf", pcaf)
        pcafile = File(pcaf)
        with open(pcaf) as jsonFile:
            pcadata = json.load(jsonFile)

    else:
        print("its in th else loop")
        response = redirect('DisplayError')
    keytotal= 0
    valuetotal=0
    for (k, v) in pcadata.items():
        print("Key: " + k)
        if(k != "PCA" and k !="warning"):
            keytotal = keytotal+int(k)
            valuetotal= valuetotal+v
            print("Value: " + str(v))

    print("keytotal",keytotal)
    print("valuetotal",valuetotal)
    finalno= round((valuetotal*100)/keytotal,2)
    finalstringno = str(valuetotal)+"/"+str(keytotal)
    print("finalno run step2: " + str(finalno))



    return JsonResponse({"foldername": foldername,"res":res,"pcajson":pcadata,"batchremovalornot":batchremovalornot,outputfile:outputfile,"finalperc":finalstringno})




def runprestep2(request,foldername):
    #print("froms tep 2 preprocessing step")
    #print(foldername)




    if foldername != '':
        print("filepath is "+foldername)

        #https://stackoverflow.com/questions/63228094/how-can-i-create-download-link-using-django
        filepath = '/Users/snehalpatil/Documents/GithubProjects/ShitingProject/tumourclass/media/'+foldername+'/cpm_input_classifier.csv'

        # Open the file for reading content
        path = open(filepath, "r")
        # Set the mime type
        mime_type, _ = mimetypes.guess_type(filepath)
        print (mime_type)
        # Set the return value of the HttpResponse
        response = HttpResponse(path, content_type=mime_type)

    #result file from the



def DisplayError(request):
    print("reached to the error")
    return render (request, "DisplayError.html")


    #return render (request, "jobsubmit.html",{"list": lines})# string of HTML code


def clustergramheatmap(request):
    print ("testng clustergram")
    folderpath="5d271d9f-966f-4252-8ed6-b9a01533932d"
    folder_save_path=os.path.join(settings.MEDIA_ROOT, "5d271d9f-966f-4252-8ed6-b9a01533932d")
    encoded_foldername = quote(folderpath)
    heatmap_json = os.path.join(settings.MEDIA_ROOT, folderpath,"mult_view.json" )
    file_save_path = os.path.join(settings.MEDIA_ROOT, folderpath,"cpm_input_classifier.csv" )

    with open(heatmap_json) as jsonFile:
        finaljson = json.load(jsonFile)


    return render (request, "clustergramheatmap.html",{"arr": file_save_path,"finaljson":finaljson,"foldername":encoded_foldername})




def runclassifier(request):

    print("****************************",request.GET.get('outputfile'))
    lines = []
    d= {}

    folderpath = request.GET.get('outputfile')
    encoded_foldername = quote(folderpath)
    batchremovalornot= request.GET.get('batchremovalornot')

    if(batchremovalornot =='combat'):
        file_save_path = os.path.join(settings.MEDIA_ROOT, folderpath,"cpm_input_classifier.csv" )
    else:
        file_save_path = os.path.join(settings.MEDIA_ROOT, folderpath,"cpm_input_classifier_no_batchre.csv" )



    #result file from the

    heatmap_file = os.path.join(settings.MEDIA_ROOT, folderpath,"heatmap_data.txt" )
    heatmap_json = os.path.join(settings.MEDIA_ROOT, folderpath,"mult_view.json" )
    pca_json = os.path.join(settings.MEDIA_ROOT, folderpath,"PCA_related.json" )
    folder_save_path=os.path.join(settings.MEDIA_ROOT, folderpath)

    print(file_save_path)

    pyfilepath = settings.PYTHONSCRIPT_PATH
    print("pyfilepath",pyfilepath)



    cols = [0, 1, 2] # add more columns here
    df = pd.DataFrame()
    models ='f,knn,gnb,svm,elasticnet'
    gene_set ='168,960'
    #def runfile(workdir, pca, input,output_dir,imumodel,gene_set):
    with open(pca_json) as pcajson:
        pcajson = json.load(pcajson)

    for key , value in pcajson.items():
        print("key is ",key)
        print("value is ",value)
        if key == "PCA":
            pcaval = value

    print(pyfilepath)
    print(pcaval)
    print(file_save_path)
    print(folder_save_path)

    runfile(pyfilepath,pcaval,file_save_path,folder_save_path,None,None)

    # result_file_path=file_save_path = os.path.join(settings.MEDIA_ROOT, folderpath,"predict_result.txt" )
    # arr = pd.read_csv(result_file_path, sep='\t', header=None, usecols=cols)
    #
    # net = Network()
    # net.load_file(heatmap_file)
    # net.cluster()
    # net.write_json_to_file('viz', heatmap_json)
    #
    #
    # with open(heatmap_json) as jsonFile:
    #     finaljson = json.load(jsonFile)


    #return redirect(request, "displayResult", {"arr": arr,"finaljson":finaljson,"foldername":encoded_foldername} )
    return redirect("displayResults", folderpath=folderpath)



def displayResults(request, folderpath):
    print("****************************",folderpath)
    encoded_foldername = quote(folderpath)
    result_file_path=os.path.join(settings.MEDIA_ROOT, folderpath,"predict_result.txt" )
    heatmap_file = os.path.join(settings.MEDIA_ROOT, folderpath,"heatmap_data.txt" )
    heatmap_json = os.path.join(settings.MEDIA_ROOT, folderpath,"mult_view.json" )
    pca_json = os.path.join(settings.MEDIA_ROOT, folderpath,"PCA_related.json" )
    folder_save_path=os.path.join(settings.MEDIA_ROOT, folderpath)


    df = pd.DataFrame()
    arr = pd.read_csv(result_file_path, sep='\t', skiprows=1)

    net = Network()
    net.load_file(heatmap_file)
    net.cluster()
    net.write_json_to_file('viz', heatmap_json)


    with open(heatmap_json) as jsonFile:
        finaljson = json.load(jsonFile)


    return render(request, "displayResults.html", {"arr": arr,"finaljson":finaljson,"foldername":encoded_foldername} )



