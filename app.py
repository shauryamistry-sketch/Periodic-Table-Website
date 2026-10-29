
#http://127.0.0.1:5000
from flask import Flask, render_template, request
from elements import get_ele
from chemistry_engine import pred_com
app = Flask(__name__)
@app.route("/")
def home():
    all_ele = get_ele()
    selnum = request.args.get("element", type=int)
    selel = next((element for element in all_ele if element["number"]==selnum), None)
    numa = request.args.get("elea", type=int)
    numb = request.args.get("eleb", type=int)
    elea=next((element for element in all_ele if element["number"] == numa), None)
    eleb=next((element for element in all_ele if element["number"] == numb), None)
    result = None
    if elea and eleb:
        result = pred_com(elea, eleb)
    forms = result["formulas"] if result else []
    bondtype = result["bond"] if result else None
    return render_template("index.html", elements=all_ele, selel=selel, elea=elea, eleb=eleb, forms=forms, bondtype=bondtype)
if __name__=="__main__":
    app.run(debug=True)

