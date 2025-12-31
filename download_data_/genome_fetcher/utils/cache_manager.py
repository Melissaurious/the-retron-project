"""
Cache management for tracking downloads and checking file completeness
"""
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime

logger = logging.getLogger(__name__)



import hashlib
import pickle
from pathlib import Path
from typing import Optional
import logging

logger = logging.getLogger(__name__)


class ResponseCache:
    """Simple disk-based cache for HTTP responses"""
    
    def __init__(self, enabled: bool = True, cache_dir: Path = None):
        """
        Initialize response cache
        
        Args:
            enabled: Whether caching is enabled
            cache_dir: Directory to store cache files (default: ~/.genome_fetcher_cache)
        """
        self.enabled = enabled
        
        if cache_dir is None:
            self.cache_dir = Path.home() / '.genome_fetcher_cache'
        else:
            self.cache_dir = Path(cache_dir)
        
        if self.enabled:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Response cache enabled: {self.cache_dir}")
    
    def _get_cache_key(self, url: str) -> str:
        """Generate cache key from URL"""
        return hashlib.md5(url.encode()).hexdigest()
    
    def get(self, url: str) -> Optional[str]:
        """Get cached response if available"""
        if not self.enabled:
            return None
        
        cache_file = self.cache_dir / f"{self._get_cache_key(url)}.cache"
        if cache_file.exists():
            try:
                with open(cache_file, 'rb') as f:
                    cached_data = pickle.load(f)
                    logger.debug(f"Cache HIT: {url[:100]}")
                    return cached_data
            except Exception as e:
                logger.warning(f"Failed to read cache: {e}")
                return None
        
        logger.debug(f"Cache MISS: {url[:100]}")
        return None
    
    def set(self, url: str, response_text: str):
        """Cache successful response"""
        if not self.enabled:
            return
        
        cache_file = self.cache_dir / f"{self._get_cache_key(url)}.cache"
        try:
            with open(cache_file, 'wb') as f:
                pickle.dump(response_text, f)
            logger.debug(f"Cached response: {url[:100]}")
        except Exception as e:
            logger.warning(f"Failed to cache response: {e}")
    
    def clear(self):
        """Clear all cached responses"""
        if not self.enabled:
            return
        
        try:
            import shutil
            shutil.rmtree(self.cache_dir)
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            logger.info("Cache cleared")
        except Exception as e:
            logger.error(f"Failed to clear cache: {e}")


class CacheManager:
    """Manage download cache and track completed downloads"""
    
    def __init__(self, output_dir: Path):
        self.output_dir = Path(output_dir)
        self.cache_file = self.output_dir / ".download_cache.json"
        self.cache = self._load_cache()
    
    def _load_cache(self) -> Dict:
        """Load cache from disk"""
        if self.cache_file.exists():
            try:
                with open(self.cache_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                logger.warning(f"Failed to load cache: {e}")
                return {}
        return {}
    
    def _save_cache(self):
        """Save cache to disk"""
        try:
            with open(self.cache_file, 'w') as f:
                json.dump(self.cache, f, indent=2)
        except Exception as e:
            logger.error(f"Failed to save cache: {e}")
    
    def is_complete(self, terminal_name: str, required_files: List[str]) -> bool:
        """
        Check if a terminal directory is complete
        
        Args:
            terminal_name: Name of the terminal directory
            required_files: List of required file names
        """
        terminal_dir = self.output_dir / terminal_name
        
        if not terminal_dir.exists():
            return False
        
        # Check if all required files exist and have content
        for filename in required_files:
            filepath = terminal_dir / filename
            if not filepath.exists():
                logger.debug(f"Missing file: {filepath}")
                return False
            if filepath.stat().st_size < 50:  # Minimum file size
                logger.debug(f"File too small: {filepath} ({filepath.stat().st_size} bytes)")
                return False
        
        logger.debug(f"Terminal {terminal_name} is complete")
        return True
    
    # def mark_complete(self, terminal_name: str, metadata: Dict):
    #     """Mark a terminal as complete in cache"""
    #     self.cache[terminal_name] = {
    #         'completed_at': datetime.now().isoformat(),
    #         'metadata': metadata
    #     }
    #     self._save_cache()
    #     logger.debug(f"Marked {terminal_name} as complete in cache")


    def mark_complete(self, terminal_name: str, metadata: Dict, 
                    files_info: Dict = None, processing_time: float = None,
                    download_source: str = None, warnings: List = None):
        """Mark a terminal as complete in cache"""
        self.cache[terminal_name] = {
            'status': 'success',
            'completed_at': datetime.now().isoformat(),
            'processing_time_seconds': processing_time,
            'metadata': metadata,
            'files_downloaded': files_info or {},
            'download_source': download_source,
            'errors': [],
            'warnings': warnings or []
        }
        self._save_cache()
    
    # def mark_failed(self, terminal_name: str, error: str):
    #     """Mark a terminal as failed in cache"""
    #     self.cache[terminal_name] = {
    #         'failed_at': datetime.now().isoformat(),
    #         'error': error
    #     }
    #     self._save_cache()
    #     logger.debug(f"Marked {terminal_name} as failed in cache")



    def mark_failed(self, terminal_name: str, error: str, metadata: Dict = None,
                    files_info: Dict = None, processing_time: float = None,
                    error_type: str = "UnknownError", step: str = None):
        """Mark a terminal as failed in cache"""
        self.cache[terminal_name] = {
            'status': 'failed',
            'completed_at': datetime.now().isoformat(),
            'processing_time_seconds': processing_time,
            'metadata': metadata or {},
            'files_downloaded': files_info or {},
            'download_source': None,
            'errors': [{
                'type': error_type,
                'message': error,
                'timestamp': datetime.now().isoformat(),
                'step': step
            }],
            'warnings': []
        }
        self._save_cache()



    def mark_partial(self, terminal_name: str, metadata: Dict, 
                    files_info: Dict, processing_time: float,
                    errors: List, warnings: List, download_source: str = None):
        """Mark a terminal as partially complete"""
        self.cache[terminal_name] = {
            'status': 'partial',
            'completed_at': datetime.now().isoformat(),
            'processing_time_seconds': processing_time,
            'metadata': metadata,
            'files_downloaded': files_info,
            'download_source': download_source,
            'errors': errors,
            'warnings': warnings
        }
        self._save_cache()


    
    def is_cached(self, terminal_name: str) -> bool:
        """Check if terminal is in cache"""
        return terminal_name in self.cache
    
    def get_cached_info(self, terminal_name: str) -> Optional[Dict]:
        """Get cached information for a terminal"""
        return self.cache.get(terminal_name)


class DownloadLogger:
    """Logger for individual terminal downloads"""
    
    def __init__(self, log_file: Path):
        self.log_file = Path(log_file)
        self.entries = []
    
    def log(self, message: str, level: str = "INFO"):
        """Add a log entry"""
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        entry = f"[{timestamp}] [{level}] {message}"
        self.entries.append(entry)
        
        # Also write to file immediately
        with open(self.log_file, 'a') as f:
            f.write(entry + '\n')
    
    def info(self, message: str):
        self.log(message, "INFO")
    
    def warning(self, message: str):
        self.log(message, "WARNING")
    
    def error(self, message: str):
        self.log(message, "ERROR")
    
    def success(self, message: str):
        self.log(message, "SUCCESS")
    
    def get_log(self) -> str:
        """Get full log as string"""
        return '\n'.join(self.entries)


def is_file_complete(filepath: Path, min_size: int = 50) -> bool:
    """Check if a file exists and has minimum content"""
    if not filepath.exists():
        return False
    return filepath.stat().st_size >= min_size