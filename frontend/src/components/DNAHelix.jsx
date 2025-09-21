const DNAHelix = () => {
  return (
    <div className="absolute left-0 top-0 h-full w-1/3 overflow-hidden opacity-20 transform -rotate-12 -translate-x-8 dna-rotate">
      <svg
        viewBox="0 0 200 600"
        className="w-full h-full transform scale-110"
        xmlns="http://www.w3.org/2000/svg"
      >
        {/* DNA Double Helix */}
        <defs>
          <pattern id="dnaPattern" x="0" y="0" width="200" height="100" patternUnits="userSpaceOnUse">
            {/* Left strand */}
            <path
              d="M50 0 Q100 25 50 50 Q0 75 50 100"
              fill="none"
              stroke="hsl(210 60% 70%)"
              strokeWidth="3"
              strokeDasharray="4,2"
            />
            {/* Right strand */}
            <path
              d="M150 0 Q100 25 150 50 Q200 75 150 100"
              fill="none"
              stroke="hsl(210 60% 70%)"
              strokeWidth="3"
              strokeDasharray="4,2"
            />
            {/* Connecting lines */}
            <line x1="50" y1="10" x2="150" y2="10" stroke="hsl(210 40% 80%)" strokeWidth="2" />
            <line x1="70" y1="25" x2="130" y2="25" stroke="hsl(210 40% 80%)" strokeWidth="2" />
            <line x1="50" y1="40" x2="150" y2="40" stroke="hsl(210 40% 80%)" strokeWidth="2" />
            <line x1="30" y1="55" x2="170" y2="55" stroke="hsl(210 40% 80%)" strokeWidth="2" />
            <line x1="50" y1="70" x2="150" y2="70" stroke="hsl(210 40% 80%)" strokeWidth="2" />
            <line x1="70" y1="85" x2="130" y2="85" stroke="hsl(210 40% 80%)" strokeWidth="2" />

            {/* Base pairs as small circles */}
            <circle cx="50" cy="10" r="3" fill="hsl(195 100% 60%)" />
            <circle cx="150" cy="10" r="3" fill="hsl(195 100% 60%)" />
            <circle cx="70" cy="25" r="3" fill="hsl(220 90% 60%)" />
            <circle cx="130" cy="25" r="3" fill="hsl(220 90% 60%)" />
            <circle cx="50" cy="40" r="3" fill="hsl(195 100% 60%)" />
            <circle cx="150" cy="40" r="3" fill="hsl(195 100% 60%)" />
            <circle cx="30" cy="55" r="3" fill="hsl(220 90% 60%)" />
            <circle cx="170" cy="55" r="3" fill="hsl(220 90% 60%)" />
            <circle cx="50" cy="70" r="3" fill="hsl(195 100% 60%)" />
            <circle cx="150" cy="70" r="3" fill="hsl(195 100% 60%)" />
            <circle cx="70" cy="85" r="3" fill="hsl(220 90% 60%)" />
            <circle cx="130" cy="85" r="3" fill="hsl(220 90% 60%)" />
          </pattern>
        </defs>

        <rect width="200" height="1200" fill="url(#dnaPattern)" />
      </svg>
    </div>
  );
};

export default DNAHelix;