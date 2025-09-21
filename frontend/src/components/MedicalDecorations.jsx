import { Pill, Circle } from "lucide-react";

const MedicalDecorations = () => {
  return (
    <div className="absolute inset-0 pointer-events-none z-10">
      {/* Right side decorations */}
      <div className="absolute top-24 right-12 animate-pulse">
        <Pill className="w-16 h-16 text-blue-500 rotate-45 opacity-40" />
      </div>

      <div className="absolute top-1/3 right-20 animate-bounce">
        <Circle className="w-20 h-20 text-teal-400 opacity-30" />
      </div>

      <div className="absolute top-1/2 right-8 animate-pulse" style={{ animationDelay: '1s' }}>
        <Pill className="w-12 h-12 text-blue-600 rotate-12 opacity-35" />
      </div>

      <div className="absolute top-2/3 right-16 animate-bounce" style={{ animationDelay: '0.5s' }}>
        <Circle className="w-14 h-14 text-teal-500 opacity-25" />
      </div>

      <div className="absolute bottom-32 right-12 animate-pulse" style={{ animationDelay: '2s' }}>
        <Pill className="w-10 h-10 text-blue-500 -rotate-30 opacity-40" />
      </div>

      <div className="absolute bottom-20 right-24 animate-bounce" style={{ animationDelay: '1.5s' }}>
        <Circle className="w-8 h-8 text-teal-400 opacity-35" />
      </div>

      {/* Small scattered circles */}
      <div className="absolute top-16 right-32">
        <div className="w-3 h-3 rounded-full bg-blue-500 opacity-50"></div>
      </div>

      <div className="absolute top-3/4 right-6">
        <div className="w-2 h-2 rounded-full bg-teal-500 opacity-60"></div>
      </div>

      <div className="absolute bottom-40 right-28">
        <div className="w-4 h-4 rounded-full bg-blue-400 opacity-40"></div>
      </div>
    </div>
  );
};

export default MedicalDecorations;