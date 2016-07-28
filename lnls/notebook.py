#!/usr/bin/env python3

from IPython.display import HTML as _HTML
from matplotlib import rcParams as _rcParams

def turn_code_on_off_html():
    string  = '''<script type="text/javascript">
code_hide=true;
function code_toggle() {
 if (code_hide){
 $('div.output_prompt').css("color","#FFFFFF");
 $('div.input').hide();
 } else {
 $('div.output_prompt').css("color","#8b0000");
 $('div.input').show();
 }
 code_hide = !code_hide
}
$( document ).ready(code_toggle)
</script>
<script>
// define a handler
function doc_keyUp(e) {
    // this would test for the 'h' key, the 'ctrl' key and the 'shift' key at the same time
    if (e.altKey && e.ctrlKey && e.keyCode == 72) {
        code_toggle();
    }
}
// register the handler
// document.removeEventListener('keyup', doc_keyUp, false);
document.addEventListener('keyup', doc_keyUp, false);
</script>
The raw code for this IPython notebook is by default hidden for easier reading.
<form action="javascript:code_toggle()"><input type="submit"
value="Click here or press CRTL+ALT+H to toggle the code on/off."></form>'''
    return _HTML(string)

def matplotlib_notebook(): _rcParams['backend'] = 'nbAgg'

def matplotlib_inline(): _rcParams['backend'] = 'module://ipykernel.pylab.backend_inline'
