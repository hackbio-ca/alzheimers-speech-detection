import Card from './Card.jsx';

function Cards() {
  return (
    <div id="cards-section">
      <p className="cards-title">Our Solution Features</p>
      <div id="cards">
        <Card text="Machine Learning"/>
        <Card text="Big Data"/>
        <Card text="Neural Networks"/>
      </div>
    </div>
  );
}

export default Cards;