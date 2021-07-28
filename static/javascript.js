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

// validation function for index
function validateIndex() {
    'use strict';
    window.addEventListener('load', function() {
        // Fetch all the forms we want to apply custom Bootstrap validation styles to
        var forms = document.getElementsByClassName('needs-validation');
        // Loop over them and prevent submission
        var validation = Array.prototype.filter.call(forms, function(form) {
            form.addEventListener('input', function(event) {
                if (form.checkValidity() === false) {
                    event.preventDefault();
                    event.stopPropagation();
                }
                form.classList.add('was-validated');
            }, false);
            form.addEventListener('submit', function(event) {
                if (form.checkValidity() === false) {
                  event.preventDefault();
                  event.stopPropagation();
                }
                form.classList.add('was-validated');
            }, false);
        });
  }, false);
}



function validateListForm() {
    var file = document.forms["varForm"]["file"];
    if (file == null) {
        alert("Please provide a file of variants");
        return false;
    }
}


function isAlphaNumeric(str) {
  var code, i, len;

  for (i = 0, len = str.length; i < len; i++) {
    code = str.charCodeAt(i);
    if (!(code > 47 && code < 58) && // numeric (0-9)
        !(code > 64 && code < 91) && // upper alpha (A-Z)
        !(code > 96 && code < 123)) { // lower alpha (a-z)
      return false;
    }
  }
  return true;
};