from flask import Flask, request, jsonify
from flask_cors import CORS

app = Flask(__name__)
CORS(app)

@app.route("/api/find", methods=["POST"])
def find_side_effects():
    # data is a string value representing the user inputted medication name
    # it should be passed to the database to make sure its valid
    # if valid, pass it to the model and retrieve some value
    # ensure that the retrieved value is in json format before returning it
    # once returned, the frontend should format the data in a way that looks good to the user
    data = request.get_json()
    query = data.get("query", "")
    return jsonify({"result": f"You searched for: {query}"})

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5001, debug=True)