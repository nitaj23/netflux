import os

from sqlalchemy import create_engine, event
from sqlalchemy.orm import sessionmaker

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

# ------------------------------------------------------------------
# URL resolution
# Default: SQLite file next to the project root.
# Override: set DATABASE_URL=postgresql+psycopg2://... for Supabase.
# ------------------------------------------------------------------
DATABASE_URL: str = os.getenv("DATABASE_URL", "sqlite:///./mfghub.db")

_is_sqlite = DATABASE_URL.startswith("sqlite")
_connect_args = {"check_same_thread": False} if _is_sqlite else {}

engine = create_engine(
    DATABASE_URL,
    connect_args=_connect_args,
    # Surface SQL to stdout when DB_ECHO=1; off by default.
    echo=os.getenv("DB_ECHO", "").lower() in ("1", "true"),
)

# ------------------------------------------------------------------
# SQLite tuning: WAL journal for better concurrent reads and
# foreign-key enforcement (SQLite disables FK checks by default).
# These pragmas are no-ops on PostgreSQL.
# ------------------------------------------------------------------
if _is_sqlite:
    @event.listens_for(engine, "connect")
    def _set_sqlite_pragmas(dbapi_conn, _connection_record):
        cur = dbapi_conn.cursor()
        cur.execute("PRAGMA journal_mode=WAL")
        cur.execute("PRAGMA foreign_keys=ON")
        cur.execute("PRAGMA synchronous=NORMAL")
        cur.close()

# ------------------------------------------------------------------
# Session factory
# expire_on_commit=False: prevents DetachedInstanceError when ORM
# objects are accessed after the session that loaded them is closed
# (common in FastAPI response serialization).
# ------------------------------------------------------------------
SessionLocal = sessionmaker(
    bind=engine,
    autocommit=False,
    autoflush=False,
    expire_on_commit=False,
)

# Alias used by scripts (export_registry.py, populate_db.py) that do:
#   with Session() as db: ...
Session = SessionLocal


def init_db() -> None:
    """Create all tables that don't exist yet. Safe to call repeatedly."""
    from src.schema import Base  # deferred to avoid circular import at module load
    Base.metadata.create_all(bind=engine)
