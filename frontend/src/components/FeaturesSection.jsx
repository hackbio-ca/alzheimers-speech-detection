const FeaturesSection = () => {
  return (
    <section className="features">
      <div className="container features__grid">
        {/* 1 */}
        <article className="feature">
          <svg className="feature__icon" viewBox="0 0 24 24" aria-hidden="true">
            <circle cx="11" cy="11" r="7" stroke="currentColor" strokeWidth="2" fill="none"/>
            <line x1="16.65" y1="16.65" x2="21" y2="21" stroke="currentColor" strokeWidth="2" />
          </svg>
          <h3 className="feature__title">Comprehensive Search</h3>
          <p className="feature__body">
            Query a wide range of medications to get a complete view of potential adverse effects.
          </p>
        </article>

        {/* 2 */}
        <article className="feature">
          <svg className="feature__icon" viewBox="0 0 24 24" aria-hidden="true">
            <path d="M12 3v6l4 2" stroke="currentColor" strokeWidth="2" fill="none" strokeLinecap="round"/>
            <circle cx="12" cy="12" r="9" stroke="currentColor" strokeWidth="2" fill="none"/>
          </svg>
          <h3 className="feature__title">AI-Driven Insights</h3>
          <p className="feature__body">
            Leveraging neural networks and large datasets, HackBio offers up-to-date, accurate results.
          </p>
        </article>

        {/* 3 */}
        <article className="feature">
          <svg className="feature__icon" viewBox="0 0 24 24" aria-hidden="true">
            <circle cx="6" cy="6" r="2" fill="currentColor"/>
            <circle cx="18" cy="6" r="2" fill="currentColor"/>
            <circle cx="12" cy="18" r="2" fill="currentColor"/>
            <path d="M8 7.5l8 0M6.8 7.9l4.6 8.2M17.2 7.9l-4.6 8.2" stroke="currentColor" strokeWidth="2" fill="none"/>
          </svg>
          <h3 className="feature__title">User-Friendly Platform</h3>
          <p className="feature__body">
            HackBio's interface is designed for ease of use, suitable for both researchers and general users.
          </p>
        </article>
      </div>
    </section>
  );
};

export default FeaturesSection;