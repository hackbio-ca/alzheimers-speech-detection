function Find() {
    const scrollToSection = (id) => {
        const element = document.getElementById(id);
        if (element) {
            const yOffset = -80;
            const y = element.getBoundingClientRect().top + window.scrollY + yOffset;
            window.scrollTo({ top: y, behavior: "smooth" });
        }
    };

    return (
    <>
    <div id="find">
        <p id="problemStatement">How will my drug affect gene expression?</p>
        <div id="findButtons">
            <button id="bodyButton" onClick={() => scrollToSection("about")}>Learn more</button>
        </div>
    </div>
    </>
    )
}

export default Find;