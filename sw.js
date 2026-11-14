// ホーム画面にアプリとして入れるためのサービスワーカー（何も保存せず、通信はそのまま通す）
// 2026-11-06: マイ圃場の水管理のお知らせ（入水したまま・中干しが長い）を、アプリを閉じていてもスマホに通知する（push）。押すと、ふくい圃場マップを開く
self.addEventListener('install',e=>{ self.skipWaiting(); });
self.addEventListener('activate',e=>{ e.waitUntil(self.clients.claim()); });
self.addEventListener('fetch',e=>{ });
// badge: スマホの上の段（通知の小さいアイコン）。Android は形（透明度）だけを使うので、色付きのアイコンでは白い■になる → マーク（3×3の角丸の四角）を白・透明地で（2026-11-07 利用者の指摘）
self.addEventListener('push',e=>{ let d={}; try{ d=e.data?e.data.json():{}; }catch(x){ d={body:e.data?e.data.text():''}; }
  e.waitUntil(self.registration.showNotification(String(d.title||'ふくい圃場マップ').slice(0,80),{body:String(d.body||'').slice(0,600),icon:'icon-192.png',badge:'badge-96.png',tag:String(d.tag||'fhm'),renotify:true,data:{url:'./'}})); });
self.addEventListener('notificationclick',e=>{ e.notification.close(); const url=new URL('./',self.registration.scope).href;
  e.waitUntil(self.clients.matchAll({type:'window',includeUncontrolled:true}).then(cs=>{ for(const c of cs){ if(c.url.startsWith(self.registration.scope)&&'focus' in c) return c.focus(); } return self.clients.openWindow(url); })); });
