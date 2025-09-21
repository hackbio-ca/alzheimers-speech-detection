import { useState } from "react";
import { Pill, Search, Loader2 } from "lucide-react";

const SearchSection = () => {
  const [searchTerm, setSearchTerm] = useState("");
  const [result, setResult] = useState("");
  const [isLoading, setIsLoading] = useState(false);

  const handleSearch = async () => {
    if (searchTerm.trim()) {
      setIsLoading(true);
      try {
        const res = await fetch("http://localhost:5001/api/find", {
          method: "POST",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ query: searchTerm }),
        });
        const data = await res.json();
        setResult(data.result);
      } catch (err) {
        console.error(err);
        setResult("Error fetching results");
      } finally {
        setIsLoading(false);
      }
    }
  };

  const handleKeyPress = (e) => {
    if (e.key === "Enter") {
      handleSearch();
    }
  };

  return (
    <main className="flex-1 flex items-center justify-center px-6 py-12 relative">
      <div className="w-fit mx-auto">
        {/* Main Card Container */}
        <div
          className="p-16 border-4 border-red-500 rounded-[20px] bg-[#FFFFFF] shadow-[0_8px_24px_rgba(10,40,80,0.12)] w-[600px] h-[600px]"
          style={{ fontFamily: "Inter, system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif" }}
        >
          {/* Title Row - MOVED UP SLIGHTLY */}
          <div className="flex items-center justify-center gap-3" style={{ marginTop: '80px' }}>
            <Pill className="w-6 h-6 text-[#0E2A47] flex-shrink-0" strokeWidth={2} />
            <h2
              className="text-[28px] leading-[34px] font-[800] text-[#0E2A47] m-0"
              style={{ letterSpacing: "-0.2px" }}
            >
              Find Drug Side Effects
            </h2>
          </div>

          {/* Description - MOVED CLOSER TO TITLE */}
          <p
            className="text-[16px] leading-[24px] text-[#2B3A4A] m-0 text-center max-w-[54ch] mx-auto"
            style={{ marginTop: '16px' }}
          >
            Discover potential side effects of common drugs, powered by AI and biomedical data.
          </p>

          {/* Search Field */}
          <div className="flex justify-center" style={{ marginTop: '32px' }}>
            <div
              className="relative flex items-center bg-[#F7FAFC] border border-[rgba(14,42,71,0.15)] rounded-[12px] px-4 h-[64px] focus-within:border-[#0AA0A8] focus-within:shadow-[0_0_0_3px_rgba(10,160,168,0.35)] w-[480px] box-border"
            >
              <input
                type="text"
                placeholder="e.g. aspirin, ibuprofen"
                value={searchTerm}
                onChange={(e) => setSearchTerm(e.target.value)}
                onKeyPress={handleKeyPress}
                aria-label="Search for drug side effects"
                className="flex-1 border-0 outline-0 bg-transparent text-[16px] leading-[22px] text-[#1E293B] placeholder-[rgba(30,41,59,0.45)] text-center pr-8"
                style={{ fontFamily: "inherit" }}
              />
              <Search className="absolute right-3 w-5 h-5 text-[rgba(43,58,74,0.7)] pointer-events-none" strokeWidth={2} />
            </div>
          </div>

          {/* SMALLER INVISIBLE SPACER */}
          <div style={{ height: '40px' }}></div>

          {/* Primary Button */}
          <div className="flex justify-center">
            <button
              onClick={handleSearch}
              disabled={isLoading}
              className="w-[300px] h-[56px] rounded-[14px] bg-gradient-to-r from-[#0E2A6B] to-[#0AA0A8] text-[#FFFFFF] font-[700] text-[18px] leading-[22px] border-0 cursor-pointer hover:brightness-[1.03] active:translate-y-[1px] focus:outline-none focus:shadow-[0_0_0_3px_rgba(10,160,168,0.35)] disabled:opacity-50 disabled:cursor-not-allowed flex items-center justify-center gap-2"
              style={{ fontFamily: "inherit" }}
            >
              {isLoading && <Loader2 className="w-5 h-5 animate-spin" />}
              {isLoading ? "Searching..." : "Search"}
            </button>
          </div>

          {/* Results */}
          {result && (
            <div className="text-center" style={{ marginTop: '32px' }}>
              <p className="text-[14px] leading-[20px] text-[rgba(43,58,74,0.65)] m-0">
                <strong className="font-[600]">Results:</strong> {result}
              </p>
            </div>
          )}

          {/* Examples */}
          <div className="text-center" style={{ marginTop: '32px' }}>
            <p className="text-[14px] leading-[20px] text-[rgba(43,58,74,0.65)] m-0">
              <strong className="font-[600]">Examples:</strong>{" "}
              <span className="font-[600] text-[#3C4C60]">Aspirin</span>{" "}
              <span className="text-[#0E2A6B]">→</span> nausea, bleeding, stomach irritation
            </p>
            <p className="text-[14px] leading-[20px] text-[rgba(43,58,74,0.65)] m-0 mt-1">
              <span className="text-[rgba(43,58,74,0.65)]">·</span>{" "}
              <span className="font-[600] text-[#3C4C60]">Ibuprofen</span>{" "}
              <span className="text-[#0E2A6B]">→</span> dizziness, rash
            </p>
          </div>
        </div>
      </div>
    </main>
  );
};

export default SearchSection;
