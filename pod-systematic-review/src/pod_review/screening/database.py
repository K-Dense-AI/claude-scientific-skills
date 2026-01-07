"""Database operations for screening workflow."""

import logging
import sqlite3
from datetime import datetime
from pathlib import Path
from typing import Any

from ..models import Citation, ScreeningDecision, ScreeningPhase, ScreeningRecord

logger = logging.getLogger(__name__)


class ScreeningDatabase:
    """Manage screening decisions in SQLite database."""

    def __init__(self, db_path: str | Path) -> None:
        """
        Initialize screening database.

        Args:
            db_path: Path to SQLite database file
        """
        self.db_path = Path(db_path)
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        self._init_db()

    def _init_db(self) -> None:
        """Initialize database schema."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Citations table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS citations (
                id TEXT PRIMARY KEY,
                title TEXT NOT NULL,
                authors TEXT,
                year INTEGER,
                journal TEXT,
                abstract TEXT,
                doi TEXT,
                pmid TEXT,
                source_database TEXT,
                data_json TEXT,
                import_date TIMESTAMP
            )
        """)

        # Screening records table
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS screening_records (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                citation_id TEXT NOT NULL,
                reviewer_id TEXT NOT NULL,
                phase TEXT NOT NULL,
                decision TEXT NOT NULL,
                exclusion_reason TEXT,
                notes TEXT,
                timestamp TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                FOREIGN KEY (citation_id) REFERENCES citations(id),
                UNIQUE(citation_id, reviewer_id, phase)
            )
        """)

        # Screening progress tracking
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS screening_progress (
                citation_id TEXT PRIMARY KEY,
                phase TEXT NOT NULL,
                status TEXT NOT NULL,
                conflict BOOLEAN DEFAULT 0,
                adjudicated BOOLEAN DEFAULT 0,
                final_decision TEXT,
                FOREIGN KEY (citation_id) REFERENCES citations(id)
            )
        """)

        conn.commit()
        conn.close()
        logger.info(f"Database initialized at {self.db_path}")

    def add_citations(self, citations: list[Citation]) -> None:
        """
        Add citations to database.

        Args:
            citations: List of citations to add
        """
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        for citation in citations:
            cursor.execute(
                """
                INSERT OR REPLACE INTO citations
                (id, title, authors, year, journal, abstract, doi, pmid, source_database, data_json, import_date)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
                (
                    citation.id,
                    citation.title,
                    "; ".join(citation.authors),
                    citation.year,
                    citation.journal,
                    citation.abstract,
                    citation.doi,
                    citation.pmid,
                    citation.source_database,
                    citation.model_dump_json(),
                    citation.import_date.isoformat(),
                ),
            )

            # Initialize progress
            cursor.execute(
                """
                INSERT OR IGNORE INTO screening_progress (citation_id, phase, status)
                VALUES (?, ?, ?)
            """,
                (citation.id, "title_abstract", "pending"),
            )

        conn.commit()
        conn.close()
        logger.info(f"Added {len(citations)} citations to database")

    def get_citation(self, citation_id: str) -> Citation | None:
        """Get citation by ID."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        cursor.execute("SELECT data_json FROM citations WHERE id = ?", (citation_id,))
        row = cursor.fetchone()
        conn.close()

        if row:
            import json

            return Citation(**json.loads(row[0]))
        return None

    def get_citations_for_screening(
        self, reviewer_id: str, phase: ScreeningPhase, limit: int = 100
    ) -> list[Citation]:
        """
        Get citations that need screening by a reviewer.

        Args:
            reviewer_id: Reviewer ID
            phase: Screening phase
            limit: Maximum number to return

        Returns:
            List of citations needing review
        """
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Get citations not yet reviewed by this reviewer in this phase
        cursor.execute(
            """
            SELECT c.data_json
            FROM citations c
            LEFT JOIN screening_records sr
                ON c.id = sr.citation_id
                AND sr.reviewer_id = ?
                AND sr.phase = ?
            WHERE sr.id IS NULL
            LIMIT ?
        """,
            (reviewer_id, phase.value, limit),
        )

        citations = []
        for row in cursor.fetchall():
            import json

            citations.append(Citation(**json.loads(row[0])))

        conn.close()
        return citations

    def add_screening_record(self, record: ScreeningRecord) -> None:
        """
        Add screening decision.

        Args:
            record: Screening record
        """
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        cursor.execute(
            """
            INSERT OR REPLACE INTO screening_records
            (citation_id, reviewer_id, phase, decision, exclusion_reason, notes, timestamp)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
            (
                record.citation_id,
                record.reviewer_id,
                record.phase.value,
                record.decision.value,
                record.exclusion_reason,
                record.notes,
                record.timestamp.isoformat(),
            ),
        )

        # Update progress
        self._update_progress(cursor, record.citation_id, record.phase)

        conn.commit()
        conn.close()

    def _update_progress(
        self, cursor: sqlite3.Cursor, citation_id: str, phase: ScreeningPhase
    ) -> None:
        """Update screening progress for a citation."""
        # Check if there's a conflict
        cursor.execute(
            """
            SELECT decision FROM screening_records
            WHERE citation_id = ? AND phase = ?
        """,
            (citation_id, phase.value),
        )

        decisions = [row[0] for row in cursor.fetchall()]

        if len(decisions) >= 2:
            # Check for conflict
            conflict = len(set(decisions)) > 1
            status = "conflict" if conflict else "reviewed"

            # If no conflict, determine final decision
            final_decision = None if conflict else decisions[0]

            cursor.execute(
                """
                INSERT OR REPLACE INTO screening_progress
                (citation_id, phase, status, conflict, final_decision)
                VALUES (?, ?, ?, ?, ?)
            """,
                (citation_id, phase.value, status, conflict, final_decision),
            )

    def get_conflicts(self, phase: ScreeningPhase) -> list[tuple[Citation, list[ScreeningRecord]]]:
        """
        Get citations with conflicting decisions.

        Args:
            phase: Screening phase

        Returns:
            List of (citation, screening_records) tuples
        """
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        cursor.execute(
            """
            SELECT citation_id FROM screening_progress
            WHERE phase = ? AND conflict = 1 AND adjudicated = 0
        """,
            (phase.value,),
        )

        conflicts = []
        for row in cursor.fetchall():
            citation_id = row[0]
            citation = self.get_citation(citation_id)
            if citation:
                records = self.get_screening_records(citation_id, phase)
                conflicts.append((citation, records))

        conn.close()
        return conflicts

    def get_screening_records(
        self, citation_id: str, phase: ScreeningPhase | None = None
    ) -> list[ScreeningRecord]:
        """Get all screening records for a citation."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        if phase:
            cursor.execute(
                """
                SELECT citation_id, reviewer_id, phase, decision, exclusion_reason, notes, timestamp
                FROM screening_records
                WHERE citation_id = ? AND phase = ?
            """,
                (citation_id, phase.value),
            )
        else:
            cursor.execute(
                """
                SELECT citation_id, reviewer_id, phase, decision, exclusion_reason, notes, timestamp
                FROM screening_records
                WHERE citation_id = ?
            """,
                (citation_id,),
            )

        records = []
        for row in cursor.fetchall():
            records.append(
                ScreeningRecord(
                    citation_id=row[0],
                    reviewer_id=row[1],
                    phase=ScreeningPhase(row[2]),
                    decision=ScreeningDecision(row[3]),
                    exclusion_reason=row[4],
                    notes=row[5],
                    timestamp=datetime.fromisoformat(row[6]),
                )
            )

        conn.close()
        return records

    def adjudicate_conflict(
        self, citation_id: str, phase: ScreeningPhase, decision: ScreeningDecision, adjudicator_id: str, notes: str = ""
    ) -> None:
        """
        Adjudicate a conflict.

        Args:
            citation_id: Citation ID
            phase: Screening phase
            decision: Final decision
            adjudicator_id: Adjudicator ID
            notes: Adjudication notes
        """
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        # Add adjudication record
        cursor.execute(
            """
            INSERT INTO screening_records
            (citation_id, reviewer_id, phase, decision, notes, timestamp)
            VALUES (?, ?, ?, ?, ?, ?)
        """,
            (
                citation_id,
                adjudicator_id,
                phase.value,
                decision.value,
                notes,
                datetime.now().isoformat(),
            ),
        )

        # Update progress
        cursor.execute(
            """
            UPDATE screening_progress
            SET adjudicated = 1, final_decision = ?
            WHERE citation_id = ? AND phase = ?
        """,
            (decision.value, citation_id, phase.value),
        )

        conn.commit()
        conn.close()

    def get_screening_stats(self) -> dict[str, Any]:
        """Get screening statistics."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()

        stats: dict[str, Any] = {}

        # Total citations
        cursor.execute("SELECT COUNT(*) FROM citations")
        stats["total_citations"] = cursor.fetchone()[0]

        # Progress by phase
        for phase in ["title_abstract", "fulltext"]:
            cursor.execute(
                """
                SELECT status, COUNT(*)
                FROM screening_progress
                WHERE phase = ?
                GROUP BY status
            """,
                (phase,),
            )
            phase_stats = dict(cursor.fetchall())
            stats[f"{phase}_progress"] = phase_stats

        conn.close()
        return stats
