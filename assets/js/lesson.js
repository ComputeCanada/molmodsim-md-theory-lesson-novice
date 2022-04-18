// Make all tables striped by default.
$("table").addClass("table table-striped");


// Handle foldable challenges and solutions (on click and at start).
$(".solution").click(function(event) {
    var trigger = $(event.target).has(".fold-unfold").length > 0
               || $(event.target).filter(".fold-unfold").length > 0;
    if (trigger) {
        $(">*:not(h2)", this).toggle(400);
        $(">h2>span.fold-unfold", this).toggleClass("glyphicon-collapse-down glyphicon-collapse-up");
        event.stopPropagation();
    }
});
$(".solution").each(function() {
    $(">*:not(h2)", this).toggle();
    var h2 = $("h2:first", this);
    h2.append("<span class='fold-unfold glyphicon glyphicon-collapse-down'></span>");
});


// Handle searches.
// Relies on document having 'meta' element with name 'search-domain'.
function google_search() {
  var query = document.getElementById("google-search").value;
  var domain = $("meta[name=search-domain]").attr("value");
  window.open("https://www.google.com/search?q=" + query + "+site:" + domain);
}

// function to shrink the life cycle bar when scrolling
$(function(){
    $('#life-cycle').data('size','big');
});

$(window).scroll(function(){
    if($(document).scrollTop() > 0)
    {
        if($('#life-cycle').data('size') == 'big')
        {
            $('#life-cycle').data('size','small');
            $('#life-cycle').stop().animate({
                padding: '5px'
            },100);
        }
    }
    else
    {
        if($('#life-cycle').data('size') == 'small')
        {
            $('#life-cycle').data('size','big');
            $('#life-cycle').stop().animate({
                padding: '15px'
            },100);
        }
    }
});


// Cookie functions for toggle-swicth state
function setCookie(cname, cvalue, exdays) {
  const d = new Date();
  d.setTime(d.getTime() + (exdays * 24 * 60 * 60 * 1000));
  let expires = "expires="+d.toUTCString();
  document.cookie = cname + "=" + cvalue + ";" + expires + ";path=/";
}

function getCookie(cname) {
  let name = cname + "=";
  let ca = document.cookie.split(';');
  for(let i = 0; i < ca.length; i++) {
    let c = ca[i];
    while (c.charAt(0) == ' ') {
      c = c.substring(1);
    }
    if (c.indexOf(name) == 0) {
      return c.substring(name.length, c.length);
    }
  }
  return "";
}

function init_instructor_mode_from_cookie() {
  let instructor_view_state = getCookie("instructor_view_state");
  if (instructor_view_state == "off") {
      $(".toggle_instructor_view > label > input").prop( "checked", false );
      $(".instructor_notes").hide();
      $(".self_study_text").show();
  } else if (instructor_view_state == "on") {
      $(".toggle_instructor_view > label > input").prop( "checked", true );
      $(".self_study_text").hide();
      $(".instructor_notes").show();
  } else {
    instructor_view_state = "off";
    if (instructor_view_state != "" && instructor_view_state != null) {
      setCookie("instructor_view_state", instructor_view_state, cookie_lifetime);
    }
  }

  if ( $(".self_study_text").length < 1 &&
       $(".instructor_view_state").length < 1 ) {
    $(".toggle_instructor_view").hide();
  }
}

// Initialize instructor_notes as hidden.
init_instructor_mode_from_cookie()

// Toogle between self_study_text and instructor_notes
$(".toggle_instructor_view").click(function(event) {
    if ($('.toggle_instructor_view > label > input').is(":checked")) {
        $(".self_study_text").hide();
        $(".instructor_notes").show();
        setCookie("instructor_view_state", "on", cookie_lifetime);
    } else {
        $(".self_study_text").show();
        $(".instructor_notes").hide();
        setCookie("instructor_view_state", "off", cookie_lifetime);
    }
});
