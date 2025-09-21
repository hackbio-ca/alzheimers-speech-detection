import { Button } from "./ui/button";
import { Dna } from "lucide-react";

const Header = ({ setPage }) => {
  return (
    <header className="sticky top-0 w-full px-6 py-4 bg-white/80 backdrop-blur-sm border-b border-gray-200 z-0">
      <nav className="flex items-center justify-between max-w-6xl mx-auto">
        <div className="flex items-center space-x-2">
          <Dna className="w-12 h-12 text-primary" />
          <span className="text-4xl font-bold text-primary">HackBio</span>
        </div>
        <div className="flex items-center space-x-8">
          <Button variant="nav" onClick={() => setPage("home")}>
            Home
          </Button>
          <Button variant="nav" onClick={() => setPage("findsideeffects")}>
            Find Side Effects
          </Button>
        </div>
      </nav>
    </header>
  );
};

export default Header;