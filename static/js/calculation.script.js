'use strict';

$(document).ready(function() {
  
  let accordion = $(".accordion__header");
  
  for (let i = 0; i < accordion.length; i++) {
    accordion[i].addEventListener("click", function () {
      this.classList.toggle("active");
      let panel = this.nextElementSibling;
      if (panel.style.maxHeight) {
        panel.style.maxHeight = null;
        panel.parentNode.classList.remove('isOpen')
      } else {
        panel.style.maxHeight = panel.scrollHeight + "px";
        panel.parentNode.classList.add('isOpen')
      }
    });
  }
})