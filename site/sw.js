// ホーム画面にアプリとして入れるための、からっぽのサービスワーカー（何も保存せず、通信はそのまま通す）
self.addEventListener('install',e=>{ self.skipWaiting(); });
self.addEventListener('activate',e=>{ e.waitUntil(self.clients.claim()); });
self.addEventListener('fetch',e=>{ });
