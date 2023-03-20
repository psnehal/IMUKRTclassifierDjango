function selectcpm(id)
{


    var chk=id.value;

    console.log("check val: "+chk);


    //var checkb = document.getElementById("logcpm");
    //console.log(checkb);




    if(chk == 0)
    {
        console.log("its in the no loop");
        document.getElementById("rangeCheck").innerHTML = '';
    }
    else
    {
        console.log("its in the yes loop");
        // document.getElementById("beff").disabled = false;

        var list =
         '<div  class="form-group row">'
        +'<div  id="test1" class="col-sm-10 col-md-9">'
        +'<table>'
            +'<tr>'
            +'<tr><td><input id="batchfile" type="file" name="batchfile" size="30" class="formObject"></td></tr>'
            +'</table></div></div>';




        document.getElementById("rangeCheck").innerHTML = list;
    }




}


function selectFFmethod(id)
{

    console.log(id.value);

    chk = id.value;

    if(chk == 'combat' )
    {
        var list =
        '<div  class="form-group row">'

        +'<div  id="test1" class="col-sm-10 col-md-9">'
            +'<label for="test1" class="">FF method</label>'
        +'<table>'
        +'<tr><td><input type="radio" name="rc1" value="yes"/><span class="formText"> Yes </span></td></tr>'
        +'<tr><td><input type="radio" name="rc2" value="no" ><span class="formText"> No </span></td></tr>'
        +'<tr><td> Upload a meta file include batch effect you want to remove.</td></tr>'
        +'<tr><td><input id="test2" type="file" name="test2meta" size="30" class="formObject"></td></tr>'
        +'</table>';







        document.getElementById("ffmethod").innerHTML = list;

    }
    else
    {
        document.getElementById("ffmethod").innerHTML = '';

    }


}



function expandAdvancedOptions()
{
    console.log('in the funciton');
    document.getElementById("more").innerHTML =
        '<span class="text">The five ML methods used are random forest, k-nearest neighbors, '
       +'Gaussian Naïve Bayes, support vector machine, and elastic net logistic regression. Each can be '
       +'performed with 3 different gene sets as features containing 7, 168, and 960 genes, respectively.'
       +'By default, all ML methods and the 2 larger gene sets are used. The 7 gene set is only recommended'
       +'if the full transcriptome is not available.'
       +'</br></br>Batch effects, if present, are adjusted using ComBat. The classifier has been validated using RNA-seq'
       +'data from a mixture of fresh frozen and FFPE tumor samples. Compared to the original clustering-based'
       +'subtype definition, it correctly classifies all original 18 UM and 66 TCGA samples tested.'
       +'</br></br>More detailed methods are available here.'
       +'</span> <a href = "#" onclick="closeExp()">less</a>';


}

function closeExp()
{
    document.getElementById("more").innerHTML = '<a href = "#" onclick="expandAdvancedOptions()">more</a>';
}
