import Home from './Home.jsx'
import FindSideEffects from './FindSideEffects.jsx'
import { useState } from "react";

function App() {
    const handleHomeClick = () => {
    if (page === "home") {
      // already on home -> scroll to top
      window.scrollTo({ top: 0, behavior: "smooth" });
    } else {
      // on another page -> go back to home
      setPage("home");
    }
  };

   const scrollToSection = (id) => {
    const element = document.getElementById(id);
    if (element) {
      const yOffset = -80; // adjust for sticky header if needed
      const y = element.getBoundingClientRect().top + window.scrollY + yOffset;
      window.scrollTo({ top: y, behavior: "smooth" });
    }
  };

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

  const [page, setPage] = useState("home");
  
  return (
  <>
  <div id="header">
    Side Effects
    <nav className="nav-links">
      <button onClick={() => handleHomeClick()}>Home</button>
      <button onClick={() => scrollToSection("about")}>About</button>
      <button onClick={() => setPage("findsideeffects")}>Find Side Effects</button>
    </nav>
  </div>
  
  {page === "home" && <Home/>}
  {page === "findsideeffects" && <FindSideEffects/>}
  </>
  );
}

export default App;