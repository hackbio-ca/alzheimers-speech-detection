import { useState } from "react";
import Footer from "./Footer";

function FindSideEffects() {
  const [input, setInput] = useState("");
  const [result, setResult] = useState("");

    const handleSubmit = async (e) => {
    e.preventDefault();
    try {
        const res = await fetch("http://localhost:5001/api/find", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ query: input }),
        });
        const data = await res.json();
        setResult(data.result);
    } catch (err) {
        console.error(err);
        setResult("Error fetching results");
    }
    };

    return (
    <>
      <div id="findSideEffects">
        <form id="sideEffectsForm" onSubmit={handleSubmit}>
          <h1>Find Effects</h1>
          <input
            type="text"
            placeholder="Enter drug candidate (e.g. Aspirin)"
            value={input}
            onChange={(e) => setInput(e.target.value)}
            required
          />
          <button id="submission" type="submit">Submit</button>
          {result && <p className="result">{result}</p>}
        </form>
      </div>
      <Footer />
    </>
  );
}

export default FindSideEffects;