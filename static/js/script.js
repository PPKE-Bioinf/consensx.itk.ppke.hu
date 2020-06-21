$(document).ready(function() {
    
  setTimeout(function () {
    $('#pageLoader').remove();
  }, 1500);
  
  let _body = $('body');
  let moodChangerButton = $('#moodChanger');
  let moodValue = localStorage.getItem('nightMode');
  if (moodValue === '0' || moodValue ===undefined) {
    _body.removeClass('nightMode');
  } else {
    _body.addClass('nightMode');
  }
  
  moodChangerButton.click(function () {
    let moodValue = localStorage.getItem('nightMode');
    if (moodValue === '1') {
      _body.removeClass('nightMode');
      $(this).html('<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="feather feather-sun"><circle cx="12" cy="12" r="5"></circle><line x1="12" y1="1" x2="12" y2="3"></line><line x1="12" y1="21" x2="12" y2="23"></line><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"></line><line x1="18.36" y1="18.36" x2="19.78" y2="19.78"></line><line x1="1" y1="12" x2="3" y2="12"></line><line x1="21" y1="12" x2="23" y2="12"></line><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"></line><line x1="18.36" y1="5.64" x2="19.78" y2="4.22"></line></svg> Night Mode OFF')
      localStorage.setItem('nightMode', '0');
    } else {
      _body.addClass('nightMode');
      $(this).html('<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="feather feather-moon"><path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"></path></svg> Night Mode ON')
      localStorage.setItem('nightMode', '1');
    }
    
  });
  
  
  let accordion = $(".accordion__header");
  
  for (let i = 0; i < accordion.length; i++) {
    accordion[i].addEventListener("click", function() {
      this.classList.toggle("active");
      let panel = this.nextElementSibling;
      if (panel.style.maxHeight){
        panel.style.maxHeight = null;
        panel.parentNode.classList.remove('isOpen')
      } else {
        panel.style.maxHeight = panel.scrollHeight + "px";
        panel.parentNode.classList.add('isOpen')
      }
    });
  }
  
  $("#pdb_file").on('change', function() {
    
    if($(this).get(0).files.length !== 0) {
      let filename = $(this).val();
      let extension = filename.replace(/^.*\./, '');
      if (extension === filename) {
        extension = '';
      } else {
        extension = extension.toLowerCase();
      }
  
      switch (extension) {
        case 'pdb':
          break;
        default:
          Swal.fire({
            title: 'Incorrect file type',
            html: '<div class="dialog__content">' +
              '<div>It seems the file you\'ve added is not a PDB file.</div></div>',
            animation: false,
            background: '#fff',
            padding: '30px 24px',
            confirmButtonText: 'Ok',
            confirmButtonColor: 'transparent',
            confirmButtonClass: 'button--dialogConfirm',
            customClass: {
              container: 'container-class',
              popup: 'dialog dialog-isVisible',
              title: 'dialog__header',
              content: 'dialog__content',
              actions: 'dialog__actions'
            },
          });
          $("#pdb_file").val('');
      }
    }
  });
  $("#xplor_file").on('change', function() {
    
    if($(this).get(0).files.length !== 0) {
      let filename = $(this).val();
      let extension = filename.replace(/^.*\./, '');
      if (extension === filename) {
        extension = '';
      } else {
        extension = extension.toLowerCase();
      }
      
      switch (extension) {
        case 'str':
          break;
        default:
          Swal.fire({
            title: 'Incorrect file type',
            html: '<div class="dialog__content">' +
              '<div>It seems the file you\'ve added is not a NMR-Star file.</div></div>',
            animation: false,
            background: '#fff',
            padding: '30px 24px',
            confirmButtonText: 'Ok',
            confirmButtonColor: 'transparent',
            confirmButtonClass: 'button--dialogConfirm',
            customClass: {
              container: 'container-class',
              popup: 'dialog dialog-isVisible',
              title: 'dialog__header',
              content: 'dialog__content',
              actions: 'dialog__actions'
            },
          });
          $("#pdb_file").val('');
      }
    }
  });
  $("#brmb_File").on('change', function() {
    
    if($(this).get(0).files.length !== 0) {
      let filename = $(this).val();
      let extension = filename.replace(/^.*\./, '');
      if (extension === filename) {
        extension = '';
      } else {
        extension = extension.toLowerCase();
      }
      
      switch (extension) {
        case 'str':
          break;
        default:
          Swal.fire({
            title: 'Incorrect file type',
            html: '<div class="dialog__content">' +
              '<div>It seems the file you\'ve added is not a NMR-Star file.</div></div>',
            animation: false,
            background: '#fff',
            padding: '30px 24px',
            confirmButtonText: 'Ok',
            confirmButtonColor: 'transparent',
            confirmButtonClass: 'button--dialogConfirm',
            customClass: {
              container: 'container-class',
              popup: 'dialog dialog-isVisible',
              title: 'dialog__header',
              content: 'dialog__content',
              actions: 'dialog__actions'
            },
          });
          $("#pdb_file").val('');
      }
    }
  });

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
  
  $('#pdb_file').on('change', function () {
    if ($(this).get(0).files.length !== 0) {
      $("label[for='pdb_file']").parent().addClass('hasContent');
      $("label[for='pdb_file']").html('<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="feather feather-check"><polyline points="20 6 9 17 4 12"></polyline></svg>'+ " <span>" +$("#bmrb_file").val() + "</span>");
    } else {
      $("label[for='pdb_file']").parent().removeClass('hasContent');
    }
  })
  $('#xplor_file').on('change', function () {
    if ($(this).get(0).files.length !== 0) {
      $("label[for='xplor_file']").parent().addClass('hasContent');
      $("label[for='xplor_file']").html('<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="feather feather-check"><polyline points="20 6 9 17 4 12"></polyline></svg>'+ " <span>" +$("#bmrb_file").val() + "</span>");
    } else {
      $("label[for='pdb_file']").parent().removeClass('hasContent');
    }
  })
  $('#bmrb_file').on('change', function () {
    if ($(this).get(0).files.length !== 0) {
      $("label[for='bmrb_file']").parent().addClass('hasContent');
      $("label[for='bmrb_file']").html('<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="feather feather-check"><polyline points="20 6 9 17 4 12"></polyline></svg>'+ " <span>" +$("#bmrb_file").val() + "</span>");
    } else {
      $("label[for='pdb_file']").parent().removeClass('hasContent');
    }
  })
  $('#bme_weights_file').on('change', function () {
    if ($(this).get(0).files.length !== 0) {
      $("label[for='bme_weights_file']").parent().addClass('hasContent');
      $("label[for='bme_weights_file']").html('<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="feather feather-check"><polyline points="20 6 9 17 4 12"></polyline></svg>'+ " <span>" +$("#bmrb_file").val() + "</span>");
    } else {
      $("label[for='pdb_file']").parent().removeClass('hasContent');
    }
  })

    // form validation
    $('#csx_form').submit(function() {
        // var hasPDB  = ($('#pdb_file').val() != "" ) || ($('#pdb_id').val() != "" );
        hasPDB  = false;
        hasMBMR = false;
        hasNOE  = false;

        if ($('#pdb_file').val()) { hasPDB = true; }
        if ($('#bmrb_file').val()) { hasMBMR = true; }
        if ($('#xplor_file').val()) { hasNOE = true; }

        console.log(hasPDB, hasMBMR, hasNOE);

        // check if necessary inputs are selected
        if (hasPDB == false || (hasMBMR == false && hasNOE == false)) {
          Swal.fire({
            title: 'Missing file and parameters',
            html: '<div class="dialog__content">' +
              '<div>You have to upload a PDB file and at least one parameter file in order to start the calculation</div></div>',
            animation: false,
            background: '#fff',
            padding: '30px 24px',
            confirmButtonText: 'Ok',
            confirmButtonColor: 'transparent',
            confirmButtonClass: 'button--dialogConfirm',
            customClass: {
              container: 'container-class',
              popup: 'dialog dialog-isVisible',
              title: 'dialog__header',
              content: 'dialog__content',
              actions: 'dialog__actions'
            },
          })
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
        $("#calculationLoader").removeClass('isHidden');
        $(".errormsg").remove();
        $("#submit_button").prop("disabled", true);
        $("#submit_button").css('cursor', 'not-allowed');
    });

    $('#submit_test').click(function() {
        // DOM and style modifications on clients side
        $(".loader").toggleClass('hidden');
        $(".loader_text").toggleClass('hidden');
        $(".errormsg").remove();
        $("#submit_button").prop("disabled", true);
        $("#submit_button").css('cursor', 'not-allowed');
    });

    $('#fit_box').change(function() {
        if (this.checked) {
            $("#fit_range").prop('disabled', false);
            $("#fit_range").focus();
        } else {
            $("#fit_range").prop('disabled', true);
          $("#fit_range").blur();
        }
    });

    $( window ).unload(function() {
        $(".loader").toggleClass('hidden');
        $(".loader_text").toggleClass('hidden');
        $(".errormsg").remove();
        $("#submit_button").prop("disabled", false);
        $("#submit_button").css('cursor', 'auto');
    });

    $('#reset_button').click(function() {
        $( '.input--file' ).each( function() {
            var $input   = $( this ),
                $label   = $input.next( 'label' );

            if( $label.attr("for") == "pdb_file" )
                $label.html( '<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="feather feather-file-plus"><path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"></path><polyline points="14 2 14 8 20 8"></polyline><line x1="12" y1="18" x2="12" y2="12"></line><line x1="9" y1="15" x2="15" y2="15"></line></svg><span>Upload PDB file</span>' );
            else if ($label.attr("for") == "xplor_file") {
              $label.html( '<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="feather feather-file-plus"><path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"></path><polyline points="14 2 14 8 20 8"></polyline><line x1="12" y1="18" x2="12" y2="12"></line><line x1="9" y1="15" x2="15" y2="15"></line></svg><span>Upload NMR-Star file</span>' );
            } else if ($label.attr("for") == "brmb_file") {
              $label.html( '<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="feather feather-file-plus"><path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"></path><polyline points="14 2 14 8 20 8"></polyline><line x1="12" y1="18" x2="12" y2="12"></line><line x1="9" y1="15" x2="15" y2="15"></line></svg><span>Upload NMR-Star file</span>' );
            } else if ($label.attr("for") == "bme_weights_file") {
              $label.html( '<svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" class="feather feather-file-plus"><path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"></path><polyline points="14 2 14 8 20 8"></polyline><line x1="12" y1="18" x2="12" y2="12"></line><line x1="9" y1="15" x2="15" y2="15"></line></svg><span>Upload BME weights file</span>' );
            }
            
        });
    });

    'use strict';

    ;( function( $, window, document, undefined )
    {
        $( '.input--file' ).each( function()
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
