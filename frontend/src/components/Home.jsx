import Find from './Find.jsx'
import Cards from './Cards.jsx'
import Footer from './Footer.jsx'
import About from './About.jsx'
import Try from './Try.jsx'

function Home({setPage}) {
    return (
    <>
    <Find/>
    <About/>
    <Cards/>
    <Try setPage={setPage}/>
    <Footer/>
    </>
    )
}

export default Home;