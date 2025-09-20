function Try({setPage}) {
    return (
    <>
    <div id="try">
        <p>Try It Out:</p>
        <button id="tryOurSolution" onClick={() => setPage("findsideeffects")}>Find Side Effects</button>
    </div>
    </>
    )
}

export default Try;