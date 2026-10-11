from flask import Flask, render_template, request
from elements import get_ele
app = Flask(__name__)
@app.route("/")
def home():
    all_ele = get_ele()
    selnum = request.args.get("element", type=int)
    selel = next((element for element in all_ele if element["number"]==selnum), None)
    return render_template("index.html", elements=all_ele, selel=selel)
if __name__=="__main__":
    app.run(debug=True)