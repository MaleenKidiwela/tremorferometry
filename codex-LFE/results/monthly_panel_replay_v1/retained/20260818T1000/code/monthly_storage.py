"""Single-writer, bounded staging ledger for monthly production integration.

This module manages file ownership and verified cleanup, not scientific validation.
It cannot turn a detector diagnostic into an approved production method.
"""
from pathlib import Path
import json,re
from .provenance import write_json,sha256

REQUIRED={'events','families','coverage','controls','provenance','verification'}

class MonthlyStorage:
    def __init__(self,root,budget_bytes):
        self.root=Path(root).resolve()
        workspace=Path(__file__).resolve().parents[1]
        if not self.root.is_relative_to(workspace):raise ValueError('Storage must be inside codex-LFE')
        if budget_bytes<=0:raise ValueError('Positive budget required')
        self.budget=int(budget_bytes)
        self.root.mkdir(parents=True,exist_ok=True)
        self.ledger=self.root/'ledger.json'
        self.state=json.loads(self.ledger.read_text()) if self.ledger.exists() else dict(active_month=None,units={})
    def save(self):write_json(self.ledger,self.state)
    def begin_month(self,month):
        if not re.fullmatch(r'\d{4}-(0[1-9]|1[0-2])',month):raise ValueError('Expected YYYY-MM')
        active=self.state['active_month']
        if active is not None and active!=month:raise ValueError('An unfinished month is already active')
        completed=self.state.get('completed_months',[])
        if month in completed:raise ValueError('Month already complete')
        if completed and month>=min(completed):raise ValueError('Months must run backward')
        self.state['active_month']=month;self.save()
    def owned(self,relative):
        p=self.root/relative
        if p.is_symlink() or not p.resolve().is_relative_to(self.root/'staged'):raise ValueError('Only owned staged files allowed')
        return p
    def begin_unit(self,unit):
        """Register a unit even when no waveform is returned.

        Empty units still require all retained products and verified coverage;
        an empty file list alone never establishes an absence of earthquakes.
        """
        if self.state['active_month'] is None:raise ValueError('No active month')
        if not re.fullmatch(r'[A-Za-z0-9_.-]+',unit):raise ValueError('Unsafe name')
        key=self.state['active_month']+'/'+unit
        record=self.state['units'].setdefault(key,dict(status='staging',files={}))
        if record['status']!='staging':raise ValueError('Unit already committed')
        self.save()
        return record
    def stage(self,unit,name,data):
        if self.state['active_month'] is None:raise ValueError('No active month')
        if not re.fullmatch(r'[A-Za-z0-9_.-]+',unit) or not re.fullmatch(r'[A-Za-z0-9_.-]+',name):raise ValueError('Unsafe name')
        key=self.state['active_month']+'/'+unit
        record=self.begin_unit(unit)
        relative=f'staged/{key}/{name}';path=self.owned(relative)
        if path.exists():raise ValueError('Staged filename already exists')
        used=sum(p.stat().st_size for p in (self.root/'staged').rglob('*') if p.is_file())
        if used+len(data)>self.budget:raise ValueError('Temporary waveform budget exceeded')
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('xb') as f:f.write(data)
        record['files'][relative]=dict(bytes=len(data),sha256=sha256(path));self.save()
        return path
    def commit(self,unit,artifacts):
        key=self.state['active_month']+'/'+unit;record=self.state['units'][key]
        if record['status']!='staging':raise ValueError('Unit must be staging')
        if not REQUIRED.issubset(artifacts):raise ValueError('Missing retained artifact categories')
        retained={}
        for kind,value in artifacts.items():
            p=Path(value).resolve()
            if not p.is_relative_to(self.root/'retained') or not p.is_file():raise ValueError('Artifacts must be retained files')
            retained[kind]=dict(path=str(p.relative_to(self.root)),sha256=sha256(p))
        verification=json.loads(Path(artifacts['verification']).read_text())
        if verification.get('passed') is not True:raise ValueError('Output verification did not pass')
        # Save an atomic commitment before any deletion is possible.
        record.update(status='committed',artifacts=retained);self.save()
    def commit_verified_unit(self,unit,artifacts,start_s,end_s,station_ids):
        """Validate experimental unit products before allowing raw cleanup."""
        from .catalog_products import write_verification
        if 'verification' not in artifacts:raise ValueError('Verification output path required')
        for value in artifacts.values():
            if not Path(value).resolve().is_relative_to(self.root/'retained'):
                raise ValueError('Artifacts must be retained files')
        report=write_verification(artifacts['verification'],artifacts,start_s,end_s,station_ids)
        self.commit(unit,artifacts)
        return report
    def cleanup(self,unit):
        key=self.state['active_month']+'/'+unit;record=self.state['units'][key]
        if record['status'] not in ['committed','cleaning']:raise ValueError('Unit is not committed')
        for item in record['artifacts'].values():
            p=self.root/item['path']
            if not p.is_file() or sha256(p)!=item['sha256']:raise ValueError('Retained output changed or missing')
        for relative,item in record['files'].items():
            p=self.owned(relative)
            if p.exists() and sha256(p)!=item['sha256']:raise ValueError('Staged data changed; preserve for inspection')
        record['status']='cleaning';self.save()
        for relative in record['files']:
            p=self.owned(relative)
            if p.exists():p.unlink()
        record['status']='cleaned';self.save()
    def finish_month(self,expected_units):
        month=self.state['active_month']
        if month is None or not expected_units:raise ValueError('Month and explicit expected units required')
        actual={key.split('/',1)[1] for key in self.state['units'] if key.startswith(month+'/')}
        if actual!=set(expected_units):raise ValueError('Work unit coverage incomplete')
        if any(self.state['units'][month+'/'+u]['status']!='cleaned' for u in expected_units):raise ValueError('Unfinished units')
        self.state.setdefault('completed_months',[]).append(month)
        self.state['active_month']=None;self.save()
