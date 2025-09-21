import Home from './Home.jsx'
import FindSideEffects from './FindSideEffects.jsx'
import { useState } from "react";

function App() {
    const handleHomeClick = () => {
    if (page === "home") {
      window.scrollTo({ top: 0, behavior: "smooth" });
    } else {
      setPage("home");
    }
  };

   const scrollToSection = (id) => {
    const element = document.getElementById(id);
    if (element) {
      const yOffset = -80;
      const y = element.getBoundingClientRect().top + window.scrollY + yOffset;
      window.scrollTo({ top: y, behavior: "smooth" });
    }
  };

  const [page, setPage] = useState("home");
  
  return (
  <>
  <div id="header">
    <b id="title">Adversis</b>
    <nav className="nav-links">
      <button onClick={() => handleHomeClick()}>Home</button>
      {page === "home" && (
        <>
          <button onClick={() => scrollToSection("about")}>About</button>
          <button class="specialButton" onClick={() => setPage("findsideeffects")}>Find Effects</button>
        </>
      )}
    </nav>
  </div>

  {page === "home" && <Home setPage={setPage} />}
  {page === "findsideeffects" && <FindSideEffects />}
  </>
  );
}

export default App;