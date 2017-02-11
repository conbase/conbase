  function change(type){ 
    var elems = document.getElementsByClassName('cell');
    for(var i=0; i < elems.length; i++){
        if (type=="min"){
            var elems2 = document.getElementsByClassName('hidden');
            for(var j=0; j < elems2.length; j++){
                
                elems2[j].style.display = 'none';
            }
            elems[i].style.width = '8em';
        } else {
            elems[i].style.width = '27em';
        }

    }
  }

  function show(id){
    change("max");
    var elems = document.getElementsByClassName(id);
    for(var i=0; i < elems.length; i++){
        elems[i].style.width = '27em';
        if (elems[i].style.display == 'block'){
                elems[i].style.display = 'none';
            }
        else {
            elems[i].style.display = 'block';
        }
    }
  }
