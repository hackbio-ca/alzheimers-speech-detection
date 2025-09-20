function App() {
  const handleClick = async () => {
    try {
      const res = await fetch("/api/hello");
      if (!res.ok) throw new Error("Network response was not ok");
      const text = await res.text();  // <- use text() for plain strings
      alert(text);
    } catch (err) {
      console.error("Failed to fetch:", err);
      alert("Error fetching message from backend");
    }
  };

  return (
    <div>
      <h1>React + Flask API Test</h1>
      <button onClick={handleClick}>Fetch Message</button>
    </div>
  );
}

export default App;