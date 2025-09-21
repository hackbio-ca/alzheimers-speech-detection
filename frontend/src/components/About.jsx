import { useState, useEffect } from "react";

const About = () => {
  const [index, setIndex] = useState(0);

  const texts = [
    "Drug side effects and efficacy depend on the interactions the chemical has with the cellular processes as much as the target site in the body.",
    "Knowledge of the extent of interactions in the cellular processes is incomplete and current pre-clinical studies miss a variety of potential interactions.",
    "Comprehensive modeling of drug-transcriptome interactions allows for a larger picture to determine potential success of treatment, negative interactions, and mitigation of risk.",
  ];

  useEffect(() => {
    const interval = setInterval(() => {
      setIndex((prevIndex) => (prevIndex + 1) % texts.length);
    }, 15000);

    return () => clearInterval(interval);
  }, [texts.length]);

  return (
    <div id="about">
      <p>{texts[index]}</p>
      <div style={{ display: "flex", justifyContent: "center", gap: "10px", marginTop: "20px" }}>
        {texts.map((_, i) => (
          <button
            key={i}
            onClick={() => setIndex(i)}
            style={{
              width: "12px",
              height: "12px",
              borderRadius: "50%",
              border: "none",
              backgroundColor: i === index ? "#004e7a" : "#d9d9d9",
              cursor: "pointer",
            }}
          ></button>
        ))}
      </div>
    </div>
  );
};

export default About;