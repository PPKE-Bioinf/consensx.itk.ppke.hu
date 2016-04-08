$(document).ready(function() {

    // parse URL variable
    function getQueryVariable(variable) {
       var query = window.location.search.substring(1);
       var vars = query.split("&");
       for (var i=0;i<vars.length;i++) {
               var pair = vars[i].split("=");
               if(pair[0] == variable){return pair[1];}
       }
       return(false);
    }

    // if SVD is enabled (url), show SVD note
    if (getQueryVariable("svd") == "1") {
        console.log("SVD enabled!");
        $(".SVDnote").removeClass("hidden");
    }

    // If comp. measure is Q-VALUE, hide non-RDC inputs
    $(".mes").on( 'change', function() {
        console.log("MEASURE changed");
        var measures = document.querySelectorAll(".options");
        var measure;

        // get comp measure
        for (i = 0; i < measures.length; i++) {
            if (measures[i].checked) {
                measure = measures[i].value;
            }
        }

        if (measure == "q-value") {
            // show Q-VALUE note
            $(".QVALnote").removeClass("hidden");

            // hide non-RDC input fields
            $(".inputrange").each(function() {
                my_id    = $( this ).attr('id');
                my_range = $( this ).attr('id').split('_');
                my_type  = my_range[0];

                if (my_type != "RDC") {
                    // clear non-RDC inputs
                    $( this ).val(0);
                    $( this ).closest(".selection-table").addClass("hidden");
                }
            })
        // if compl. measure is other than Q-VALUE
        } else {
            // hide Q-VALUE note
            $(".QVALnote").addClass("hidden");

            // show non-RDC input fields
            $(".inputrange").each(function() {
                my_id    = $( this ).attr('id');
                my_range = $( this ).attr('id').split('_');
                my_type  = my_range[0];

                if (my_type != "RDC") {
                    $( this ).closest(".selection-table").removeClass("hidden");
                }
            })
        }
    });


    $(".inputrange").on( 'change', function() {
        if (getQueryVariable("svd") == "1") {
            my_range = $( this ).attr('id').split('_');
            sel_type = my_range[0];
            sel_num  = my_range[1];
            sel_val  = $( this ).val();

            if (sel_type == "RDC") {
                $(".inputrange").each(function() {
                    my_id    = $( this ).attr('id');
                    my_range = $( this ).attr('id').split('_');
                    my_type  = my_range[0];
                    my_num   = my_range[1];

                    if (my_type == "RDC" && my_num == sel_num) {
                        $( this ).val(sel_val);
                        console.log('output[for="' + my_range + '"]');
                        $('output[for="' + my_id + '"]').prop('value', sel_val);
                    }
                })
            }
        }
    });

    $("#selection").click(function () {
        $(".involve-result").toggleClass("involve", 1000, "easeOutSine")
    });

    $('#startsel').click(function() {
        var ranges = document.querySelectorAll(".inputrange");
        var measures = document.querySelectorAll(".options");
        var min_size = $("#min_ens_size").val().trim();
        var max_size = $("#max_ens_size").val().trim();
        var overdrive = $("#overdrive").val().trim();
        var measure;
        var range_selected = false;
        var command = "";
        var i;

        // add selected measure to command
        for (i = 0; i < measures.length; i++) {
            if (measures[i].checked) {
                console.log("selected measure: " + measures[i].value);
                command += "MEASURE " + measures[i].value + "/";
                measure = measures[i].value;
            }
        }

        if (min_size != "") {
            command += "MIN_SIZE" + " " + min_size + "/";
            console.log("minimum size is: " + min_size);
        }

        if (max_size != "") {
            command += "MAX_SIZE" + " " + max_size + "/";
            console.log("minimum size is: " + max_size);
        }

        if (overdrive != "") {
            command += "OVERDRIVE" + " " + overdrive + "/";
            console.log("overdrive is: " + overdrive);
        }

        // check if measure is selected
        if (measure == undefined) {
            alert("Please select compliance measure!");
            return false;
        }

        // add involved parameters to command
        for (i = 0; i < ranges.length; i++) {
            if (ranges[i].value != 0) {
                command += ranges[i].id + " " + ranges[i].value + "/";
                range_selected = true;
            }
        }

        // check if at least one parameter is involved in selection
        if (!range_selected) {
            alert("Please include at least one parameter!");
            return false;
        }

        // UPLAD HERE!
        console.log(command);

        $('#sel_str').val(command);
    });
});
