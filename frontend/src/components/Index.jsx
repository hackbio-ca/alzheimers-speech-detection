import Header from "./Header";
import SearchSection from "./SearchSection";
import Footer from "./Footer";
import DNAHelix from "./DNAHelix";
import MedicalDecorations from "./MedicalDecorations";

const Index = ({ setPage }) => {
  return (
    <div className="min-h-screen flex flex-col subtle-gradient relative overflow-hidden">
      <DNAHelix />
      <Header setPage={setPage} />
      <MedicalDecorations />
      <main className="flex-1 flex items-center justify-center px-6 py-12 relative">
        <SearchSection />
      </main>
      <Footer setPage={setPage} />
    </div>
  );
};

export default Index;