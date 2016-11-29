$(document).ready(function() {

  var $window = $(window)
  var items = $("#localtoc a")
  if (items.length) {
    items.prepend('<i class="icon-chevron-right pull-right"></i>')

    var sections = $(".section");
    var section_id = sections[0].id;
    if (section_id.indexOf('.') >= 0) {
      items[0].href = "#" + $(sections[0]).find('>:first-child').attr('id');
    } else {
      items[0].href = "#" + section_id;
    }

    items.first().addClass("navfirst")
    items.last().addClass("navlast")
    $("#localtoc ul").addClass("nav nav-list")
    $("#localtoc > ul").addClass("toc")
  }

  setTimeout(function() {
      $('#sidebar,.subnav').affix({
        offset: {
          top: 140,
          bottom: 20
        }
      });
    }, 100);

  $('.collapse').not('#collapseToolbox').collapse('hide');

  // Alter logo based on URL
  var url = document.URL;
  var logo = document.getElementById("logo");
  if (url.search("nmwiz") > -1) {
    logo.src = "http://prody.csb.pitt.edu/_static/nmwiz.png";
  } else if (url.search("evol") > -1) {
    logo.src = "http://prody.csb.pitt.edu/_static/evol.png";
  } else if (url.search("drugui") > -1) {
    logo.src = "http://prody.csb.pitt.edu/_static/drugui_logo.png";
  }
  else if (url.search("comd") > -1) {
    logo.src = "http://prody.csb.pitt.edu/_static/comdlogo.png";
  }
  else if (url.search("memanm") > -1) {
    logo.src = "http://prody.csb.pitt.edu/_static/membranm.png";
  }
  else if (url.search("mechstiff") > -1) {
    logo.src = "http://prody.csb.pitt.edu/_static/mechstifflogo.png";
  }else {
    logo.src = "http://prody.csb.pitt.edu/_static/logo.png";
  }
  
  if ($('#homepagenav').length) {
    if (url.search("evol") > -1) {
      $('.nav-pills > li:nth-child(2)').addClass('active');
    } else if (url.search("nmwiz") > -1) {
      $('.nav-pills > li:nth-child(3)').addClass('active');
    } else if (url.search("memanm") > -1) {
      $('.nav-pills > li:nth-child(4)').addClass('active');
    } else if (url.search("mechstiff") > -1) {
      $('.nav-pills > li:nth-child(5)').addClass('active');
    } else if (url.search("drugui") > -1) {
      $('.nav-pills > li:nth-child(6)').addClass('active');
    } else if (url.search("comd") > -1) {
      $('.nav-pills > li:nth-child(7)').addClass('active');
    } else if (url.search("downloads") > -1) {
      $('.nav-pills > li:nth-child(8)').addClass('active');
    } else if (url.search("tutorials") > -1) {
      $('.nav-pills > li:nth-child(9)').addClass('active');
    } else if (url.search("statistics") > -1) {
      $('.nav-pills > li:nth-child(10)').addClass('active');
    } else {
      $('.nav-pills > li:nth-child(1)').addClass('active');
    }
  }

});
