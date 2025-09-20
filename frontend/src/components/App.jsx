import Header from './Header.jsx'
import Find from './Find.jsx'
import Cards from './Cards.jsx'
import Footer from './Footer.jsx'
import About from './About.jsx'

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
    <>
    <Header/>
    <Find/>
    <About/>
    <Cards/>
    <Footer/>
    </>
  );
}

export default App;