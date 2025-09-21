function Try({setPage}) {
    return (
    <>
    <div id="try">
        <p><b>Try It Out:</b></p>
        <button id="tryOurSolution" onClick={() => setPage("findsideeffects")}>Find Effects</button>
    </div>
    </>
    )
}

export default Try;