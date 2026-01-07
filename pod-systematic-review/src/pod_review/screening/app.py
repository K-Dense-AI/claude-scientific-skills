"""Streamlit app for citation screening."""

import streamlit as st
from pathlib import Path

from ..config import ReviewConfig, get_settings
from ..models import ScreeningDecision, ScreeningPhase, ScreeningRecord
from .database import ScreeningDatabase


def run_screening_app(config_path: str = "config.yaml") -> None:
    """
    Run the Streamlit screening application.

    Args:
        config_path: Path to configuration file
    """
    st.set_page_config(page_title="POD Systematic Review - Screening", layout="wide")

    # Load configuration
    config = ReviewConfig(config_path)
    settings = get_settings()

    # Initialize database
    db = ScreeningDatabase(settings.screening_db_path)

    # Sidebar navigation
    st.sidebar.title("Screening App")
    page = st.sidebar.radio(
        "Navigation", ["Reviewer Login", "Screen Citations", "Adjudication", "Statistics"]
    )

    # Session state for reviewer
    if "reviewer_id" not in st.session_state:
        st.session_state.reviewer_id = None
    if "current_index" not in st.session_state:
        st.session_state.current_index = 0

    if page == "Reviewer Login":
        show_login_page(config)
    elif page == "Screen Citations":
        if st.session_state.reviewer_id:
            show_screening_page(db, config)
        else:
            st.warning("Please login first")
    elif page == "Adjudication":
        if st.session_state.reviewer_id:
            show_adjudication_page(db, config)
        else:
            st.warning("Please login first")
    elif page == "Statistics":
        show_statistics_page(db)


def show_login_page(config: ReviewConfig) -> None:
    """Show reviewer login page."""
    st.title("Reviewer Login")

    reviewers = config.reviewers
    reviewer_options = [f"{r['id']} - {r['name']}" for r in reviewers]

    selected = st.selectbox("Select Reviewer", reviewer_options)

    if st.button("Login"):
        if selected:
            reviewer_id = selected.split(" - ")[0]
            st.session_state.reviewer_id = reviewer_id
            st.success(f"Logged in as {selected}")
            st.rerun()


def show_screening_page(db: ScreeningDatabase, config: ReviewConfig) -> None:
    """Show citation screening page."""
    st.title("Citation Screening")

    reviewer_id = st.session_state.reviewer_id

    # Select phase
    phase_option = st.radio("Screening Phase", ["Title/Abstract", "Full Text"])
    phase = ScreeningPhase.TITLE_ABSTRACT if phase_option == "Title/Abstract" else ScreeningPhase.FULLTEXT

    # Get citations to screen
    citations = db.get_citations_for_screening(reviewer_id, phase, limit=1000)

    if not citations:
        st.info("No citations available for screening in this phase.")
        return

    # Navigation
    col1, col2, col3 = st.columns([1, 3, 1])
    with col1:
        if st.button("← Previous") and st.session_state.current_index > 0:
            st.session_state.current_index -= 1
            st.rerun()
    with col2:
        st.write(
            f"Citation {st.session_state.current_index + 1} of {len(citations)}"
        )
    with col3:
        if st.button("Next →") and st.session_state.current_index < len(citations) - 1:
            st.session_state.current_index += 1
            st.rerun()

    # Get current citation
    if st.session_state.current_index >= len(citations):
        st.session_state.current_index = 0

    citation = citations[st.session_state.current_index]

    # Display citation
    st.subheader(citation.title)

    col1, col2 = st.columns([2, 1])

    with col1:
        st.write(f"**Authors:** {citation.get_author_string()}")
        st.write(f"**Journal:** {citation.journal} ({citation.year})")
        if citation.doi:
            st.write(f"**DOI:** {citation.doi}")
        if citation.pmid:
            st.write(f"**PMID:** {citation.pmid}")

    with col2:
        st.write(f"**Source:** {citation.source_database}")
        st.write(f"**ID:** {citation.id}")

    st.write("**Abstract:**")
    st.write(citation.abstract if citation.abstract else "_No abstract available_")

    # Screening decision
    st.markdown("---")
    st.subheader("Screening Decision")

    col1, col2 = st.columns(2)

    with col1:
        decision_option = st.radio(
            "Decision",
            ["Include", "Exclude", "Unclear"],
            key=f"decision_{citation.id}",
        )
        decision_map = {
            "Include": ScreeningDecision.INCLUDE,
            "Exclude": ScreeningDecision.EXCLUDE,
            "Unclear": ScreeningDecision.UNCLEAR,
        }
        decision = decision_map[decision_option]

    with col2:
        exclusion_reason = ""
        if decision == ScreeningDecision.EXCLUDE:
            exclusion_reasons = config.exclusion_reasons
            reason_options = [f"{r['code']}: {r['reason']}" for r in exclusion_reasons]
            selected_reason = st.selectbox("Exclusion Reason", reason_options)
            exclusion_reason = selected_reason.split(":")[0] if selected_reason else ""

    notes = st.text_area("Notes (optional)", key=f"notes_{citation.id}")

    if st.button("Submit Decision", type="primary"):
        record = ScreeningRecord(
            citation_id=citation.id,
            reviewer_id=reviewer_id,
            phase=phase,
            decision=decision,
            exclusion_reason=exclusion_reason,
            notes=notes,
        )
        db.add_screening_record(record)
        st.success("Decision recorded!")

        # Move to next citation
        if st.session_state.current_index < len(citations) - 1:
            st.session_state.current_index += 1
        st.rerun()


def show_adjudication_page(db: ScreeningDatabase, config: ReviewConfig) -> None:
    """Show adjudication page for conflicts."""
    st.title("Conflict Adjudication")

    reviewer_id = st.session_state.reviewer_id

    # Check if user is adjudicator
    adjudicators = [r["id"] for r in config.reviewers if r.get("role") == "adjudicator"]
    if reviewer_id not in adjudicators:
        st.warning("You are not authorized to adjudicate conflicts.")
        return

    # Select phase
    phase_option = st.radio("Phase", ["Title/Abstract", "Full Text"], key="adj_phase")
    phase = ScreeningPhase.TITLE_ABSTRACT if phase_option == "Title/Abstract" else ScreeningPhase.FULLTEXT

    # Get conflicts
    conflicts = db.get_conflicts(phase)

    if not conflicts:
        st.info("No conflicts to adjudicate in this phase.")
        return

    st.write(f"Total conflicts: {len(conflicts)}")

    # Show conflicts
    for idx, (citation, records) in enumerate(conflicts):
        with st.expander(f"Conflict {idx + 1}: {citation.title[:100]}..."):
            st.subheader(citation.title)
            st.write(f"**Authors:** {citation.get_author_string()}")
            st.write(f"**Journal:** {citation.journal} ({citation.year})")
            st.write("**Abstract:**")
            st.write(citation.abstract if citation.abstract else "_No abstract available_")

            st.markdown("---")
            st.write("**Reviewer Decisions:**")

            for record in records:
                st.write(
                    f"- {record.reviewer_id}: **{record.decision.value}** "
                    f"({record.exclusion_reason or 'no reason'})"
                )
                if record.notes:
                    st.write(f"  Notes: _{record.notes}_")

            # Adjudication decision
            adj_col1, adj_col2 = st.columns(2)

            with adj_col1:
                adj_decision = st.radio(
                    "Final Decision",
                    ["Include", "Exclude", "Unclear"],
                    key=f"adj_decision_{citation.id}",
                )
                decision_map = {
                    "Include": ScreeningDecision.INCLUDE,
                    "Exclude": ScreeningDecision.EXCLUDE,
                    "Unclear": ScreeningDecision.UNCLEAR,
                }
                final_decision = decision_map[adj_decision]

            with adj_col2:
                adj_notes = st.text_area("Adjudication Notes", key=f"adj_notes_{citation.id}")

            if st.button("Submit Adjudication", key=f"submit_{citation.id}"):
                db.adjudicate_conflict(citation.id, phase, final_decision, reviewer_id, adj_notes)
                st.success("Adjudication recorded!")
                st.rerun()


def show_statistics_page(db: ScreeningDatabase) -> None:
    """Show screening statistics."""
    st.title("Screening Statistics")

    stats = db.get_screening_stats()

    st.metric("Total Citations", stats.get("total_citations", 0))

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Title/Abstract Screening")
        ta_stats = stats.get("title_abstract_progress", {})
        for status, count in ta_stats.items():
            st.write(f"- {status}: {count}")

    with col2:
        st.subheader("Full-Text Screening")
        ft_stats = stats.get("fulltext_progress", {})
        if ft_stats:
            for status, count in ft_stats.items():
                st.write(f"- {status}: {count}")
        else:
            st.write("Not started")


if __name__ == "__main__":
    run_screening_app()
