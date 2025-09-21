import { useState, useEffect } from 'react';

const HeroSection = ({ setPage }) => {
  const [displayText, setDisplayText] = useState('');
  const fullText = 'Discover drug side effects powered by AI and biomedical data';

  useEffect(() => {
    let index = 0;
    const timer = setInterval(() => {
      if (index < fullText.length) {
        setDisplayText(fullText.slice(0, index + 1));
        index++;
      } else {
        clearInterval(timer);
      }
    }, 50); // Adjust speed here

    return () => clearInterval(timer);
  }, []);

  return (
    <section className="hero">
      <div className="hero__inner">
        <h1 className="hero__title">Welcome to HackBio</h1>
        <p className="hero__tagline">
          {displayText}<span className="cursor">|</span>
        </p>
        <button className="btn btn--primary" onClick={() => setPage("findsideeffects")}>Find Side Effects</button>
      </div>
    </section>
  );
};

export default HeroSection;