$(document).ready(function() {

    console.log("DOC READY!");

    console.log("PDB FILE: " + $("#pdb_file").val());

    if ($('#pdb_file').val()) {
        svg = '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="17" viewBox="0 0 20 17"><path d="M10 0l-5.2 4.9h3.3v5.1h3.8v-5.1h3.3l-5.2-4.9zm9.3 11.5l-3.2-2.1h-2l3.4 2.6h-3.5c-.1 0-.2.1-.2.1l-.8 2.3h-6l-.8-2.2c-.1-.1-.1-.2-.2-.2h-3.6l3.4-2.6h-2l-3.2 2.1c-.4.3-.7 1-.6 1.5l.6 3.1c.1.5.7.9 1.2.9h16.3c.6 0 1.1-.4 1.3-.9l.6-3.1c.1-.5-.2-1.2-.7-1.5z"/></svg>';
        $("label[for='pdb_file']").html(svg + " <span>" +$("#pdb_file").val() + "</span>");
    }

    if ($('#xplor_file').val()) {
        svg = '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="17" viewBox="0 0 20 17"><path d="M10 0l-5.2 4.9h3.3v5.1h3.8v-5.1h3.3l-5.2-4.9zm9.3 11.5l-3.2-2.1h-2l3.4 2.6h-3.5c-.1 0-.2.1-.2.1l-.8 2.3h-6l-.8-2.2c-.1-.1-.1-.2-.2-.2h-3.6l3.4-2.6h-2l-3.2 2.1c-.4.3-.7 1-.6 1.5l.6 3.1c.1.5.7.9 1.2.9h16.3c.6 0 1.1-.4 1.3-.9l.6-3.1c.1-.5-.2-1.2-.7-1.5z"/></svg>';
        $("label[for='xplor_file']").html(svg + " <span>" +$("#xplor_file").val() + "</span>");
    }

    if ($('#bmrb_file').val()) {
        svg = '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="17" viewBox="0 0 20 17"><path d="M10 0l-5.2 4.9h3.3v5.1h3.8v-5.1h3.3l-5.2-4.9zm9.3 11.5l-3.2-2.1h-2l3.4 2.6h-3.5c-.1 0-.2.1-.2.1l-.8 2.3h-6l-.8-2.2c-.1-.1-.1-.2-.2-.2h-3.6l3.4-2.6h-2l-3.2 2.1c-.4.3-.7 1-.6 1.5l.6 3.1c.1.5.7.9 1.2.9h16.3c.6 0 1.1-.4 1.3-.9l.6-3.1c.1-.5-.2-1.2-.7-1.5z"/></svg>';
        $("label[for='bmrb_file']").html(svg + " <span>" +$("#bmrb_file").val() + "</span>");
    }


    $("#fit_range").prop('disabled', true);

    // form validation
    $('#csx_form').submit(function() {
        hasPDB  = false;
        hasMBMR = false;
        hasNOE  = false;
        hasSAXS = false;

        if ($('#pdb_file').val()) { hasPDB = true; }
        if ($('#bmrb_file').val()) { hasMBMR = true; }
        if ($('#xplor_file').val()) { hasNOE = true; }
        if ($('#saxs_dat_file').val()) { hasSAXS = true; }

        console.log(hasPDB, hasMBMR, hasNOE);

        // check if necessary inputs are selected
        if (hasPDB == false || (hasMBMR == false && hasNOE == false && hasSAXS == false)) {
            alert("You have to upload a PDB file and at least one parameter file to start calculation");
            return false;
        }

        // check PDB ID input
        // if ($('#pdb_id').val() != "") {
        //     var PDB_ID   = $('#pdb_id').val()
        //     var alphanum = /((^[0-9]+[a-z]+)|(^[a-z]+[0-9]+))+[0-9a-z]+$/i;


        //     if (!PDB_ID.match(alphanum) || PDB_ID.length != 4) {
        //         alert("Invalid PDB ID");
        //         return false;
        //     }
        // }

        // check if both PDB inputs are selected
        // var hasBothPDB = ($('#pdb_file').val() != "" ) && ($('#pdb_id').val() != "" );

        // if (hasBothPDB) {
        //     alert("Please deselect one of the PDB inputs");
        //     return false;
        // }

        // check selected fit range
        var fit_range = $('#fit_range').val()

        if (fit_range != "") {
            var borders = fit_range.split("-");

            if (borders[0] >= borders[1] || borders[0] == undefined || borders[1] == undefined) {
                alert("Invalid fit range");
                return false;
            }
        }

        // DOM and style modifications on clients side
        $("#collapseOne").collapse('hide');
        $("#collapseTwo").collapse('hide');
        $("#collapseThird").collapse('hide');
        $(".loader").toggleClass('hidden');
        $(".loader_text").toggleClass('hidden');
        $(".errormsg").remove();
        $("#submit_button").prop("disabled", true);
        $("#submit_button").css('cursor', 'not-allowed');
    });

    $('#submit_test').click(function() {
        // DOM and style modifications on clients side
        $("#collapseOne").collapse('hide');
        $("#collapseTwo").collapse('hide');
        $("#collapseThird").collapse('hide');
        $(".loader").toggleClass('hidden');
        $(".loader_text").toggleClass('hidden');
        $(".errormsg").remove();
        $("#submit_button").prop("disabled", true);
        $("#submit_button").css('cursor', 'not-allowed');
    });

    $('#fit_box').change(function() {
        if (this.checked) {
            $("#fit_range").prop('disabled', false);
        } else {
            $("#fit_range").prop('disabled', true);
        }
    });

    $( window ).unload(function() {
        $("#collapseOne").collapse('hide');
        $("#collapseTwo").collapse('hide');
        $("#collapseThird").collapse('hide');
        $(".loader").toggleClass('hidden');
        $(".loader_text").toggleClass('hidden');
        $(".errormsg").remove();
        $("#submit_button").prop("disabled", false);
        $("#submit_button").css('cursor', 'auto');
    });

    $('#reset_button').click(function() {
        $( '.inputfile' ).each( function() {
            var $input   = $( this ),
                $label   = $input.next( 'label' );

            if( $label.attr("for") == "pdb_file" )
                $label.html( '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="17" viewBox="0 0 20 17"><path d="M10 0l-5.2 4.9h3.3v5.1h3.8v-5.1h3.3l-5.2-4.9zm9.3 11.5l-3.2-2.1h-2l3.4 2.6h-3.5c-.1 0-.2.1-.2.1l-.8 2.3h-6l-.8-2.2c-.1-.1-.1-.2-.2-.2h-3.6l3.4-2.6h-2l-3.2 2.1c-.4.3-.7 1-.6 1.5l.6 3.1c.1.5.7.9 1.2.9h16.3c.6 0 1.1-.4 1.3-.9l.6-3.1c.1-.5-.2-1.2-.7-1.5z"/></svg> <span>Upload PDB file&hellip;</span>' );
            else
                $label.html( '<svg xmlns="http://www.w3.org/2000/svg" width="20" height="17" viewBox="0 0 20 17"><path d="M10 0l-5.2 4.9h3.3v5.1h3.8v-5.1h3.3l-5.2-4.9zm9.3 11.5l-3.2-2.1h-2l3.4 2.6h-3.5c-.1 0-.2.1-.2.1l-.8 2.3h-6l-.8-2.2c-.1-.1-.1-.2-.2-.2h-3.6l3.4-2.6h-2l-3.2 2.1c-.4.3-.7 1-.6 1.5l.6 3.1c.1.5.7.9 1.2.9h16.3c.6 0 1.1-.4 1.3-.9l.6-3.1c.1-.5-.2-1.2-.7-1.5z"/></svg> <span>Upload BMRB file&hellip;</span>' );
        });
    });

    'use strict';

    ;( function( $, window, document, undefined )
    {
        $( '.inputfile' ).each( function()
        {
            var $input   = $( this ),
                $label   = $input.next( 'label' ),
                labelVal = $label.html();

            $input.on( 'change', function( e )
            {
                var fileName = '';

                if( this.files && this.files.length > 1 )
                    fileName = ( this.getAttribute( 'data-multiple-caption' ) ||
                                 '' ).replace( '{count}', this.files.length );
                else if( e.target.value )
                    fileName = e.target.value.split( '\\' ).pop();

                if( fileName )
                    $label.find( 'span' ).html( fileName );
                else
                    $label.html( labelVal );
            });

            // Firefox bug fix
            $input
            .on( 'focus', function(){ $input.addClass( 'has-focus' ); })
            .on( 'blur', function(){ $input.removeClass( 'has-focus' ); });
        });
    })( jQuery, window, document );

});
