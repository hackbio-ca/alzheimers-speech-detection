import Card from './Card.jsx';

function Cards() {
  return (
    <div id="cards-section">
      <p className="cards-title"><b>Our Solution Features</b></p>
      <div id="cards">
        <Card title ="Comprehensive Search" text="Query a wide range of medications to get a complete view of potential adverse effects."/>
        <Card title ="AI-Driven Insights" text="Leveraging neural networks and large datasets, Adversis offers up-to-date, accurate results."/>
        <Card title ="User-Friendly Platform" text="Adversis's interface is designed for ease of use, suitable for both researchers and general users."/>
      </div>
    </div>
  );
}

export default Cards;