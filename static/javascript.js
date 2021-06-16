function validateForm() {
    var gene = document.forms["varForm"]["gene"].value;
    if (gene == "") {
        alert("Please provide the gene name");
        return false;
    }
    var loc = document.forms["varForm"]["position"].value;
    if (loc == "") {
        alert("Please provide the cDNA coordinates");
        return false;
    }
    //var wt = document.forms["varForm"]["wt"].value;
    //var mut = document.forms["varForm"]["mut"].value;
    //if (wt == mut) {
    //    alert("The wild type and mutant alleles must be different bases");
    //    return false;
    //}

}

function validateListForm() {
    var file = document.forms["varForm"]["file"];
    if (file == null) {
        alert("Please provide a file of variants");
        return false;
    }
}